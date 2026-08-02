from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from typing import Any

from django.db.models import OuterRef, Subquery

from database.models import (
    AxisContextPresence,
    AxisFinalEvidence,
    AxisResultKind,
)

from .import_normalization import (
    AxisImportValidationError,
    build_axis_key,
    build_axis_signature,
    derive_canonical_axis_type,
    normalize_text,
)


DEFAULT_MATCH_QUERY_CHUNK_SIZE = 1000

MATCH_SCOPE_AXIS_FINAL = "axis_final"
MATCH_SCOPE_SPONGE = "sponge"
MATCH_SCOPE_BOTH = "both"
MATCH_SCOPE_ANY = "any"

ALLOWED_MATCH_SCOPES = {
    MATCH_SCOPE_AXIS_FINAL,
    MATCH_SCOPE_SPONGE,
    MATCH_SCOPE_BOTH,
    MATCH_SCOPE_ANY,
}


class RecurrentAxisProjectMatchError(ValueError):
    pass


# Compatibility alias for older imports.
AxisFinalProjectMatchError = (
    RecurrentAxisProjectMatchError
)


def _validate_positive_integer(
    value,
    *,
    field_name: str,
) -> int:
    try:
        value = int(value)
    except (TypeError, ValueError) as exc:
        raise RecurrentAxisProjectMatchError(
            f"{field_name} must be a positive integer."
        ) from exc

    if value <= 0:
        raise RecurrentAxisProjectMatchError(
            f"{field_name} must be a positive integer."
        )

    return value


def normalize_match_scope(
    value: str | None,
) -> str:
    scope = normalize_text(
        value or MATCH_SCOPE_AXIS_FINAL
    ).lower()

    if scope not in ALLOWED_MATCH_SCOPES:
        raise RecurrentAxisProjectMatchError(
            "Unsupported project match scope. "
            "Allowed values: "
            + ", ".join(
                sorted(ALLOWED_MATCH_SCOPES)
            )
            + "."
        )

    return scope


def iter_chunks(
    values: Iterable[Any],
    *,
    chunk_size: int,
):
    chunk_size = _validate_positive_integer(
        chunk_size,
        field_name="chunk_size",
    )

    chunk = []

    for value in values:
        chunk.append(value)

        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []

    if chunk:
        yield chunk


def build_workflow_axis_identity(
    record: dict,
    *,
    row_index: int,
) -> dict:
    """
    Build the same structural identity used by the recurrent Axis importer.

    Canonical signature:
        miRNA|mRNA|lncRNA|circRNA

    ``axis_type`` from the workflow CSV is not used as a matching key.
    """
    miRNA = normalize_text(
        record.get("miRNA")
    )
    mRNA = normalize_text(
        record.get("mRNA")
    )
    lncRNA = normalize_text(
        record.get("lncRNA")
    )
    circRNA = normalize_text(
        record.get("circRNA")
    )

    canonical_axis_type = (
        derive_canonical_axis_type(
            miRNA=miRNA,
            mRNA=mRNA,
            lncRNA=lncRNA,
            circRNA=circRNA,
            row_index=row_index,
        )
    )

    axis_signature = build_axis_signature(
        miRNA=miRNA,
        mRNA=mRNA,
        lncRNA=lncRNA,
        circRNA=circRNA,
    )

    return {
        "axis_signature": axis_signature,
        "axis_key": build_axis_key(
            axis_signature
        ),
        "canonical_axis_type":
            canonical_axis_type,
    }


def add_workflow_axis_identity_to_records(
    records: list[dict] | None,
) -> list[dict]:
    """
    Copy records and append CanonicalAxis identity fields.

    One invalid row does not abort matching for every other workflow row.
    """
    enriched_records = []

    for row_index, source_record in enumerate(
        records or []
    ):
        record = dict(source_record or {})

        try:
            identity = build_workflow_axis_identity(
                record,
                row_index=row_index,
            )
        except AxisImportValidationError as exc:
            enriched_records.append({
                **record,
                "axis_signature": "",
                "axis_key": "",
                "canonical_axis_type": "",
                "axis_project_match_error":
                    str(exc),
            })
            continue

        enriched_records.append({
            **record,
            **identity,
            "axis_project_match_error": None,
        })

    return enriched_records


def get_unique_axis_keys_from_records(
    records: list[dict],
) -> list[str]:
    return sorted({
        normalize_text(
            record.get("axis_key")
        )
        for record in records
        if normalize_text(
            record.get("axis_key")
        )
    })


def apply_match_scope(
    queryset,
    *,
    match_scope: str,
):
    if match_scope == MATCH_SCOPE_AXIS_FINAL:
        return queryset.filter(
            has_axis_final=True,
        )

    if match_scope == MATCH_SCOPE_SPONGE:
        return queryset.filter(
            has_sponge=True,
        )

    if match_scope == MATCH_SCOPE_BOTH:
        return queryset.filter(
            has_axis_final=True,
            has_sponge=True,
        )

    return queryset


def build_project_key(
    *,
    dataset_source: str,
    module: str,
    dataset_name: str,
    group_type: str,
    group_by: str,
) -> str:
    return "|".join([
        normalize_text(dataset_source),
        normalize_text(module),
        normalize_text(dataset_name),
        normalize_text(group_type),
        normalize_text(group_by),
    ])


def serialize_context_project_match(
    presence: AxisContextPresence,
) -> dict:
    """
    Serialize one active Axis occurrence under one dataset/grouping context.

    In the new schema, a former "project match" is represented by:
        AxisContextPresence + AxisDatasetContext.
    """
    axis = presence.axis
    context = presence.context

    dataset_name = normalize_text(
        context.dataset_metadata_id
    )
    axis_regulation = normalize_text(
        getattr(
            presence,
            "matched_axis_regulation",
            "",
        )
    )

    project_key = build_project_key(
        dataset_source=context.dataset_source,
        module=context.module,
        dataset_name=dataset_name,
        group_type=context.group_type,
        group_by=context.group_by,
    )

    return {
        "project_key": project_key,

        "canonical_axis_id": axis.id,
        "axis_key": axis.axis_key,
        "axis_signature": axis.axis_signature,
        "axis_type": axis.axis_type,

        "context_presence_id": presence.id,
        "context_id": context.id,

        "dataset_source":
            context.dataset_source,
        "module": context.module,
        "dataset_name": dataset_name,
        "group_type": context.group_type,
        "group_by": context.group_by,

        "annotation_dir_name":
            context.annotation_dir_name,
        "annotation_file_prefix":
            context.annotation_file_prefix,

        "observation_count":
            presence.observation_count,
        "axis_final_observation_count":
            presence.axis_final_observation_count,
        "sponge_observation_count":
            presence.sponge_observation_count,

        "has_axis_final":
            presence.has_axis_final,
        "has_sponge":
            presence.has_sponge,
        "has_both_results": bool(
            presence.has_axis_final
            and presence.has_sponge
        ),

        "axis_regulation":
            axis_regulation,

        # Transitional aliases for existing frontend code.
        "project_id": context.id,
        "occurrence_id": presence.id,
        "source": context.dataset_source,
        "axis_id": axis.id,
    }


def query_project_matches_by_axis_keys(
    *,
    axis_keys: list[str],
    match_scope: str = MATCH_SCOPE_AXIS_FINAL,
    chunk_size: int = DEFAULT_MATCH_QUERY_CHUNK_SIZE,
) -> dict[str, list[dict]]:
    """
    Query active context-level project matches by CanonicalAxis.axis_key.
    """
    match_scope = normalize_match_scope(
        match_scope
    )
    chunk_size = _validate_positive_integer(
        chunk_size,
        field_name="chunk_size",
    )

    unique_axis_keys = sorted({
        normalize_text(axis_key)
        for axis_key in axis_keys
        if normalize_text(axis_key)
    })

    if not unique_axis_keys:
        return {}

    matches_by_axis_key = defaultdict(list)

    axis_regulation_subquery = (
        AxisFinalEvidence.objects
        .filter(
            observation__axis_id=OuterRef(
                "axis_id"
            ),
            observation__artifact__context_id=(
                OuterRef("context_id")
            ),
            observation__artifact__is_active=True,
            observation__artifact__context__is_active=True,
            observation__artifact__result_kind=(
                AxisResultKind.AXIS_FINAL
            ),
        )
        .order_by("observation_id")
        .values("axis_regulation")[:1]
    )

    for axis_key_chunk in iter_chunks(
        unique_axis_keys,
        chunk_size=chunk_size,
    ):
        queryset = (
            AxisContextPresence.objects
            .filter(
                axis__axis_key__in=(
                    axis_key_chunk
                ),
                context__is_active=True,
            )
            .select_related(
                "axis",
                "context",
                "context__dataset_metadata",
            )
            .annotate(
                matched_axis_regulation=Subquery(
                    axis_regulation_subquery
                )
            )
        )

        queryset = apply_match_scope(
            queryset,
            match_scope=match_scope,
        )

        queryset = queryset.order_by(
            "axis__axis_key",
            "context__dataset_source",
            "context__module",
            "context__dataset_metadata_id",
            "context__group_type",
            "context__group_by",
            "context_id",
        )

        for presence in queryset.iterator(
            chunk_size=chunk_size
        ):
            matches_by_axis_key[
                presence.axis.axis_key
            ].append(
                serialize_context_project_match(
                    presence
                )
            )

    return dict(matches_by_axis_key)


def build_axis_project_match_summary(
    *,
    records: list[dict],
    matches_by_axis_key: dict[
        str,
        list[dict],
    ],
    match_scope: str,
) -> dict:
    matched_record_count = 0
    invalid_record_count = 0

    valid_axis_keys = set()
    matched_axis_keys = set()

    matched_project_keys = set()
    matched_dataset_keys = set()
    axis_project_pairs = set()

    project_axis_keys = defaultdict(set)

    for record in records:
        axis_key = normalize_text(
            record.get("axis_key")
        )

        if not axis_key:
            invalid_record_count += 1
            continue

        valid_axis_keys.add(axis_key)

        matches = matches_by_axis_key.get(
            axis_key,
            [],
        )

        if matches:
            matched_record_count += 1
            matched_axis_keys.add(axis_key)

        for match in matches:
            project_key = match["project_key"]

            matched_project_keys.add(
                project_key
            )
            matched_dataset_keys.add((
                match["dataset_source"],
                match["dataset_name"],
            ))
            axis_project_pairs.add((
                axis_key,
                project_key,
            ))
            project_axis_keys[
                project_key
            ].add(axis_key)

    project_hits = []

    for (
        project_key,
        project_matched_axis_keys,
    ) in project_axis_keys.items():
        first_match = next(
            match
            for matches in matches_by_axis_key.values()
            for match in matches
            if match["project_key"] == project_key
        )

        project_hits.append({
            "project_key": project_key,
            "context_id":
                first_match["context_id"],
            "dataset_source":
                first_match["dataset_source"],
            "source":
                first_match["dataset_source"],
            "module":
                first_match["module"],
            "dataset_name":
                first_match["dataset_name"],
            "group_type":
                first_match["group_type"],
            "group_by":
                first_match["group_by"],
            "matched_axis_count": len(
                project_matched_axis_keys
            ),
        })

    project_hits.sort(
        key=lambda item: (
            -item["matched_axis_count"],
            item["dataset_source"],
            item["module"],
            item["dataset_name"],
            item["group_type"],
            item["group_by"],
        )
    )

    total_record_count = len(records)
    total_unique_axis_count = len(
        valid_axis_keys
    )
    matched_unique_axis_count = len(
        matched_axis_keys
    )

    return {
        "match_scope": match_scope,

        "total_record_count":
            total_record_count,
        "invalid_record_count":
            invalid_record_count,

        "total_unique_axis_count":
            total_unique_axis_count,
        "matched_record_count":
            matched_record_count,
        "unmatched_record_count": (
            total_record_count
            - matched_record_count
        ),
        "matched_unique_axis_count":
            matched_unique_axis_count,
        "unmatched_unique_axis_count": (
            total_unique_axis_count
            - matched_unique_axis_count
        ),

        "matched_dataset_count": len(
            matched_dataset_keys
        ),
        "matched_context_count": len(
            matched_project_keys
        ),
        "axis_context_match_count": len(
            axis_project_pairs
        ),
        "context_hits": project_hits,

        # Transitional aliases.
        "total_axis_count":
            total_record_count,
        "matched_axis_count":
            matched_record_count,
        "unmatched_axis_count": (
            total_record_count
            - matched_record_count
        ),
        "matched_project_count": len(
            matched_project_keys
        ),
        "project_hits": project_hits,
    }


def attach_project_matches_to_axis_records(
    *,
    records: list[dict] | None,
    max_matches_per_axis: int | None = None,
    match_scope: str = MATCH_SCOPE_AXIS_FINAL,
    chunk_size: int = DEFAULT_MATCH_QUERY_CHUNK_SIZE,
) -> dict:
    """
    Attach recurrent-reference project/context matches to workflow rows.

    Workflow rows are read-only and are not imported into recurrent tables.
    """
    match_scope = normalize_match_scope(
        match_scope
    )

    if max_matches_per_axis is not None:
        max_matches_per_axis = (
            _validate_positive_integer(
                max_matches_per_axis,
                field_name=(
                    "max_matches_per_axis"
                ),
            )
        )

    records_with_identity = (
        add_workflow_axis_identity_to_records(
            records
        )
    )

    axis_keys = get_unique_axis_keys_from_records(
        records_with_identity
    )

    matches_by_axis_key = (
        query_project_matches_by_axis_keys(
            axis_keys=axis_keys,
            match_scope=match_scope,
            chunk_size=chunk_size,
        )
    )

    summary = build_axis_project_match_summary(
        records=records_with_identity,
        matches_by_axis_key=matches_by_axis_key,
        match_scope=match_scope,
    )

    enriched_records = []

    for record in records_with_identity:
        axis_key = normalize_text(
            record.get("axis_key")
        )

        matches = matches_by_axis_key.get(
            axis_key,
            [],
        )

        full_match_count = len(matches)

        if max_matches_per_axis is None:
            response_matches = matches
        else:
            response_matches = matches[
                :max_matches_per_axis
            ]

        dataset_match_count = len({
            (
                match["dataset_source"],
                match["dataset_name"],
            )
            for match in matches
        })

        canonical_axis_id = (
            matches[0]["canonical_axis_id"]
            if matches
            else None
        )

        enriched_records.append({
            **record,

            "canonical_axis_id":
                canonical_axis_id,

            "reference_dataset_match_count":
                dataset_match_count,
            "reference_context_match_count":
                full_match_count,
            "reference_context_matches":
                response_matches,
            "reference_matches_truncated": (
                len(response_matches)
                < full_match_count
            ),

            # Transitional aliases.
            "dataset_project_match_count":
                full_match_count,
            "dataset_project_matches":
                response_matches,
        })

    return {
        "records": enriched_records,
        "summary": summary,
    }


def enrich_axis_final_response_with_project_matches(
    *,
    response_data: dict,
    max_matches_per_axis: int | None = None,
    match_scope: str = MATCH_SCOPE_AXIS_FINAL,
    chunk_size: int = DEFAULT_MATCH_QUERY_CHUNK_SIZE,
) -> dict:
    """
    Public service entry point for workflow Axis Final response enrichment.
    """
    if not isinstance(response_data, dict):
        raise RecurrentAxisProjectMatchError(
            "response_data must be a dict."
        )

    records = response_data.get(
        "results",
        [],
    )

    if not isinstance(records, list):
        raise RecurrentAxisProjectMatchError(
            "response_data['results'] "
            "must be a list."
        )

    match_scope = normalize_match_scope(
        match_scope
    )

    match_result = (
        attach_project_matches_to_axis_records(
            records=records,
            max_matches_per_axis=(
                max_matches_per_axis
            ),
            match_scope=match_scope,
            chunk_size=chunk_size,
        )
    )

    return {
        **response_data,
        "results": match_result["records"],

        "axis_reference_match_summary":
            match_result["summary"],
        "axis_reference_match_enabled":
            True,

        # Transitional aliases.
        "axis_project_match_summary":
            match_result["summary"],
        "axis_project_match_enabled":
            True,
    }
