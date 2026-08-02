from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping

from django.db import IntegrityError, transaction

from database.models import (
    AxisDatasetContext,
    AxisFinalEvidence,
    AxisObservation,
    AxisResultArtifact,
    AxisResultKind,
    CanonicalAxis,
    DatasetMetadata,
    SpongeEvidence,
)

from .import_adapters import get_axis_result_adapter
from .import_context import validate_context_spec
from .import_contracts import (
    AxisContextSpec,
    NormalizedAxisRow,
    ParsedAxisArtifact,
)
from .import_normalization import AxisImportValidationError
from .import_presence import rebuild_axis_context_presence


CANONICAL_AXIS_SIGNATURE_VERSION = 2

AXIS_FINAL_EVIDENCE_FIELDS = (
    "axis_regulation",
    "mRNA_log2FC",
    "mRNA_regulation",
    "miRNA_log2FC",
    "miRNA_regulation",
    "lncRNA_log2FC",
    "lncRNA_regulation",
    "circRNA_log2FC",
    "circRNA_regulation",
)

SPONGE_EVIDENCE_FIELDS = (
    "cor",
    "pcor",
    "mscor",
)


def import_axis_result_file(
    *,
    context_spec: AxisContextSpec,
    result_kind: str,
    file_path: str | Path,
    schema_version: str = "v1",
    dry_run: bool = False,
    skip_unchanged: bool = True,
    batch_size: int = 1000,
) -> dict:
    """
    Validate and import one Axis result artifact.

    Validation and CSV parsing happen outside the database transaction. The
    database transaction only contains context/artifact switching and row
    persistence, which keeps lock duration bounded.

    One active artifact is retained per:
        (context, result_kind)

    Re-import behavior:
        - unchanged file: return skipped=True;
        - changed file: deactivate the previous active artifact and insert a
          new artifact with observations/evidence;
        - failure: the transaction rolls back and the old artifact remains
          active.
    """
    _validate_batch_size(batch_size)
    validate_context_spec(context_spec)

    dataset_metadata = _get_dataset_metadata(
        dataset_name=context_spec.dataset_name,
    )

    adapter = get_axis_result_adapter(result_kind)

    parsed = adapter.parse_file(
        file_path=Path(file_path),
        schema_version=schema_version,
    )

    if parsed.result_kind != result_kind:
        raise AxisImportValidationError(
            "Adapter result kind does not match the requested result kind: "
            f"requested={result_kind!r}, parsed={parsed.result_kind!r}."
        )

    if dry_run:
        return {
            "success": True,
            "dry_run": True,
            "skipped": False,
            "dataset_source": context_spec.dataset_source,
            "module": context_spec.module,
            "dataset_name": context_spec.dataset_name,
            "group_type": context_spec.group_type,
            "group_by": context_spec.group_by,
            "result_kind": parsed.result_kind,
            "schema_version": parsed.schema_version,
            "file_name": parsed.file_path.name,
            "file_path": str(parsed.file_path),
            "file_sha256": parsed.file_sha256,
            "row_count": parsed.row_count,
        }

    return _persist_axis_result_file(
        context_spec=context_spec,
        dataset_metadata=dataset_metadata,
        parsed=parsed,
        skip_unchanged=skip_unchanged,
        batch_size=batch_size,
    )


def _get_dataset_metadata(
    *,
    dataset_name: str,
) -> DatasetMetadata:
    try:
        return DatasetMetadata.objects.get(
            dataset=dataset_name,
        )
    except DatasetMetadata.DoesNotExist as exc:
        raise AxisImportValidationError(
            "DatasetMetadata does not contain the dataset required by the "
            f"Axis context: {dataset_name!r}."
        ) from exc
    except DatasetMetadata.MultipleObjectsReturned as exc:
        # dataset is expected to be unique, so this indicates schema/data
        # corruption rather than a normal validation failure.
        raise RuntimeError(
            "Multiple DatasetMetadata rows were found for "
            f"dataset={dataset_name!r}."
        ) from exc


@transaction.atomic
def _persist_axis_result_file(
    *,
    context_spec: AxisContextSpec,
    dataset_metadata: DatasetMetadata,
    parsed: ParsedAxisArtifact,
    skip_unchanged: bool,
    batch_size: int,
) -> dict:
    context, context_created, context_updated_fields = upsert_context(
        context_spec=context_spec,
        dataset_metadata=dataset_metadata,
    )

    # Locking the context serializes imports for both result kinds under the
    # same context. This also closes the "no active artifact row exists yet"
    # race that cannot be protected by locking only AxisResultArtifact rows.
    context = (
        AxisDatasetContext.objects
        .select_for_update()
        .get(pk=context.pk)
    )

    active_artifacts = list(
        AxisResultArtifact.objects
        .select_for_update()
        .filter(
            context=context,
            result_kind=parsed.result_kind,
            is_active=True,
        )
        .order_by("id")
    )

    unchanged_artifact = _find_unchanged_artifact(
        active_artifacts=active_artifacts,
        parsed=parsed,
        skip_unchanged=skip_unchanged,
    )

    if unchanged_artifact is not None:
        artifact_updated_fields = _refresh_unchanged_artifact_metadata(
            artifact=unchanged_artifact,
            parsed=parsed,
        )

        return {
            "success": True,
            "dry_run": False,
            "skipped": True,
            "reason": "unchanged_file_sha256",
            "dataset_source": context.dataset_source,
            "module": context.module,
            "dataset_name": context.dataset_name,
            "group_type": context.group_type,
            "group_by": context.group_by,
            "context_id": context.id,
            "context_created": context_created,
            "context_updated_fields": context_updated_fields,
            "artifact_id": unchanged_artifact.id,
            "artifact_updated_fields": artifact_updated_fields,
            "result_kind": parsed.result_kind,
            "schema_version": parsed.schema_version,
            "file_name": parsed.file_path.name,
            "file_path": str(parsed.file_path),
            "file_sha256": parsed.file_sha256,
            "row_count": parsed.row_count,
        }

    replaced_artifact_ids = [
        artifact.id
        for artifact in active_artifacts
    ]

    if replaced_artifact_ids:
        AxisResultArtifact.objects.filter(
            id__in=replaced_artifact_ids,
        ).update(
            is_active=False,
        )

    artifact = AxisResultArtifact.objects.create(
        context=context,
        result_kind=parsed.result_kind,
        file_name=parsed.file_path.name,
        file_path=str(parsed.file_path),
        file_sha256=parsed.file_sha256,
        row_count=parsed.row_count,
        schema_version=parsed.schema_version,
        is_active=True,
    )

    axes_by_key = ensure_canonical_axes(
        rows=parsed.rows,
        batch_size=batch_size,
    )

    observations_by_row = create_observations(
        artifact=artifact,
        rows=parsed.rows,
        axes_by_key=axes_by_key,
        batch_size=batch_size,
    )

    evidence_count = create_evidence(
        result_kind=parsed.result_kind,
        rows=parsed.rows,
        observations_by_row=observations_by_row,
        batch_size=batch_size,
    )

    presence_result = rebuild_axis_context_presence(
        context=context,
        batch_size=max(batch_size, 1000),
    )

    return {
        "success": True,
        "dry_run": False,
        "skipped": False,
        "dataset_source": context.dataset_source,
        "module": context.module,
        "dataset_name": context.dataset_name,
        "group_type": context.group_type,
        "group_by": context.group_by,
        "context_id": context.id,
        "context_created": context_created,
        "context_updated_fields": context_updated_fields,
        "artifact_id": artifact.id,
        "replaced_artifact_ids": replaced_artifact_ids,
        "result_kind": parsed.result_kind,
        "schema_version": parsed.schema_version,
        "file_name": parsed.file_path.name,
        "file_path": str(parsed.file_path),
        "file_sha256": parsed.file_sha256,
        "row_count": parsed.row_count,
        "canonical_axis_count": len(axes_by_key),
        "observation_count": len(observations_by_row),
        "evidence_count": evidence_count,
        "presence": presence_result,
    }


def upsert_context(
    *,
    context_spec: AxisContextSpec,
    dataset_metadata: DatasetMetadata,
) -> tuple[AxisDatasetContext, bool, list[str]]:
    """
    Get or create the immutable context identity and refresh mutable metadata.

    Context identity:
        dataset_source
        module
        dataset_metadata
        group_type
        group_by

    Mutable fields:
        annotation_dir_name
        annotation_file_prefix
        is_active
    """
    lookup = {
        "dataset_source": context_spec.dataset_source,
        "module": context_spec.module,
        "dataset_metadata": dataset_metadata,
        "group_type": context_spec.group_type,
        "group_by": context_spec.group_by,
    }

    context = (
        AxisDatasetContext.objects
        .select_for_update()
        .filter(**lookup)
        .first()
    )

    created = False

    if context is None:
        # The nested atomic block provides a savepoint. If a concurrent worker
        # creates the same unique context first, only this create is rolled
        # back and the outer import transaction remains usable.
        try:
            with transaction.atomic():
                context = AxisDatasetContext.objects.create(
                    **lookup,
                    annotation_dir_name=(
                        context_spec.annotation_dir_name
                    ),
                    annotation_file_prefix=(
                        context_spec.annotation_file_prefix
                    ),
                    is_active=True,
                )
                created = True
        except IntegrityError:
            context = (
                AxisDatasetContext.objects
                .select_for_update()
                .get(**lookup)
            )

    updated_fields: list[str] = []

    # Empty values in a manually constructed AxisContextSpec are treated as
    # "not supplied" so they do not erase discovery metadata already stored.
    if (
        context_spec.annotation_dir_name
        and context.annotation_dir_name
        != context_spec.annotation_dir_name
    ):
        context.annotation_dir_name = (
            context_spec.annotation_dir_name
        )
        updated_fields.append("annotation_dir_name")

    if (
        context_spec.annotation_file_prefix
        and context.annotation_file_prefix
        != context_spec.annotation_file_prefix
    ):
        context.annotation_file_prefix = (
            context_spec.annotation_file_prefix
        )
        updated_fields.append("annotation_file_prefix")

    if not context.is_active:
        context.is_active = True
        updated_fields.append("is_active")

    if updated_fields:
        context.save(
            update_fields=[
                *updated_fields,
                "updated_at",
            ]
        )

    return context, created, updated_fields


def _find_unchanged_artifact(
    *,
    active_artifacts: list[AxisResultArtifact],
    parsed: ParsedAxisArtifact,
    skip_unchanged: bool,
) -> AxisResultArtifact | None:
    if not skip_unchanged:
        return None

    if len(active_artifacts) != 1:
        return None

    artifact = active_artifacts[0]

    if artifact.file_sha256 != parsed.file_sha256:
        return None

    if artifact.schema_version != parsed.schema_version:
        return None

    if artifact.row_count != parsed.row_count:
        return None

    # A matching hash with a mismatched observation count indicates an old
    # incomplete/corrupt import and must be repaired rather than skipped.
    if artifact.observations.count() != parsed.row_count:
        return None

    if parsed.result_kind == AxisResultKind.AXIS_FINAL:
        evidence_count = AxisFinalEvidence.objects.filter(
            observation__artifact=artifact,
        ).count()

    elif parsed.result_kind == AxisResultKind.SPONGE:
        evidence_count = SpongeEvidence.objects.filter(
            observation__artifact=artifact,
        ).count()

    if evidence_count != parsed.row_count:
        return None

    return artifact


def _refresh_unchanged_artifact_metadata(
    *,
    artifact: AxisResultArtifact,
    parsed: ParsedAxisArtifact,
) -> list[str]:
    updated_fields: list[str] = []
    desired_path = str(parsed.file_path)

    if artifact.file_name != parsed.file_path.name:
        artifact.file_name = parsed.file_path.name
        updated_fields.append("file_name")

    if artifact.file_path != desired_path:
        artifact.file_path = desired_path
        updated_fields.append("file_path")

    if updated_fields:
        artifact.save(
            update_fields=[
                *updated_fields,
                "updated_at",
            ]
        )

    return updated_fields


def ensure_canonical_axes(
    *,
    rows: Iterable[NormalizedAxisRow],
    batch_size: int,
) -> dict[str, CanonicalAxis]:
    rows_by_key = {
        row.axis_key: row
        for row in rows
    }

    if not rows_by_key:
        return {}

    axis_keys = list(rows_by_key)

    existing = CanonicalAxis.objects.in_bulk(
        axis_keys,
        field_name="axis_key",
    )

    missing = [
        CanonicalAxis(
            signature_version=(
                CANONICAL_AXIS_SIGNATURE_VERSION
            ),
            axis_key=row.axis_key,
            axis_signature=row.axis_signature,
            axis_type=row.axis_type,
            miRNA=row.miRNA,
            mRNA=row.mRNA,
            lncRNA=row.lncRNA,
            circRNA=row.circRNA,
        )
        for axis_key, row in rows_by_key.items()
        if axis_key not in existing
    ]

    if missing:
        CanonicalAxis.objects.bulk_create(
            missing,
            batch_size=batch_size,
            ignore_conflicts=True,
        )

    axes = CanonicalAxis.objects.in_bulk(
        axis_keys,
        field_name="axis_key",
    )

    missing_after_insert = (
        set(axis_keys) - set(axes)
    )

    if missing_after_insert:
        raise RuntimeError(
            "CanonicalAxis rows are missing after bulk creation: "
            f"{sorted(missing_after_insert)[:5]!r}."
        )

    for axis_key, axis in axes.items():
        row = rows_by_key[axis_key]

        actual = (
            axis.signature_version,
            axis.axis_signature,
            axis.axis_type,
            axis.miRNA,
            axis.mRNA,
            axis.lncRNA,
            axis.circRNA,
        )
        expected = (
            CANONICAL_AXIS_SIGNATURE_VERSION,
            row.axis_signature,
            row.axis_type,
            row.miRNA,
            row.mRNA,
            row.lncRNA,
            row.circRNA,
        )

        if actual != expected:
            raise RuntimeError(
                "Axis key collision or inconsistent canonical "
                f"normalization for axis_key={axis_key!r}; "
                f"stored={actual!r}, incoming={expected!r}."
            )

    return axes


def create_observations(
    *,
    artifact: AxisResultArtifact,
    rows: Iterable[NormalizedAxisRow],
    axes_by_key: Mapping[str, CanonicalAxis],
    batch_size: int,
) -> dict[int, AxisObservation]:
    rows = tuple(rows)

    if not rows:
        return {}

    observations = []

    for row in rows:
        try:
            axis = axes_by_key[row.axis_key]
        except KeyError as exc:
            raise RuntimeError(
                "CanonicalAxis was not resolved for "
                f"row_index={row.row_index}, "
                f"axis_key={row.axis_key!r}."
            ) from exc

        observations.append(
            AxisObservation(
                artifact=artifact,
                axis=axis,
                row_index=row.row_index,
                source_axis_id=row.source_axis_id,
                source_axis_type=row.source_axis_type,
                extra_data=dict(row.extra_data),
            )
        )

    AxisObservation.objects.bulk_create(
        observations,
        batch_size=batch_size,
    )

    # Query back by artifact rather than relying on backend-specific PK
    # population behavior of bulk_create().
    persisted = list(
        AxisObservation.objects
        .filter(artifact=artifact)
        .only("id", "row_index")
        .order_by("row_index")
    )

    observations_by_row = {
        observation.row_index: observation
        for observation in persisted
    }

    expected_row_indexes = {
        row.row_index
        for row in rows
    }

    if set(observations_by_row) != expected_row_indexes:
        raise RuntimeError(
            "Persisted AxisObservation rows do not match parsed rows; "
            f"expected={len(expected_row_indexes)}, "
            f"persisted={len(observations_by_row)}."
        )

    return observations_by_row


def create_evidence(
    *,
    result_kind: str,
    rows: Iterable[NormalizedAxisRow],
    observations_by_row: Mapping[int, AxisObservation],
    batch_size: int,
) -> int:
    rows = tuple(rows)

    if not rows:
        return 0

    if result_kind == AxisResultKind.AXIS_FINAL:
        objects = [
            AxisFinalEvidence(
                observation=observations_by_row[
                    row.row_index
                ],
                **_select_evidence_fields(
                    evidence=row.evidence,
                    allowed_fields=(
                        AXIS_FINAL_EVIDENCE_FIELDS
                    ),
                    row_index=row.row_index,
                ),
            )
            for row in rows
        ]

        AxisFinalEvidence.objects.bulk_create(
            objects,
            batch_size=batch_size,
        )

        return len(objects)

    if result_kind == AxisResultKind.SPONGE:
        objects = [
            SpongeEvidence(
                observation=observations_by_row[
                    row.row_index
                ],
                **_select_evidence_fields(
                    evidence=row.evidence,
                    allowed_fields=(
                        SPONGE_EVIDENCE_FIELDS
                    ),
                    row_index=row.row_index,
                ),
            )
            for row in rows
        ]

        SpongeEvidence.objects.bulk_create(
            objects,
            batch_size=batch_size,
        )

        return len(objects)

    raise AxisImportValidationError(
        f"Unsupported Axis result kind: {result_kind!r}."
    )


def _select_evidence_fields(
    *,
    evidence: Mapping,
    allowed_fields: tuple[str, ...],
    row_index: int,
) -> dict:
    unexpected_fields = (
        set(evidence) - set(allowed_fields)
    )

    if unexpected_fields:
        raise AxisImportValidationError(
            f"Row {row_index}: unexpected evidence fields: "
            f"{sorted(unexpected_fields)!r}."
        )

    return {
        field_name: evidence.get(field_name)
        for field_name in allowed_fields
    }


def _validate_batch_size(batch_size: int) -> None:
    if not isinstance(batch_size, int):
        raise TypeError("batch_size must be an integer.")

    if batch_size <= 0:
        raise ValueError("batch_size must be greater than zero.")
