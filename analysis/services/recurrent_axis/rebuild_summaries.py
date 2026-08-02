from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from typing import Any

from django.db import transaction
from django.db.models import Count, F, Q, QuerySet

from database.models import (
    AxisContextPresence,
    AxisDatasetSource,
    AxisFinalEvidence,
    AxisFinalRecurrentSummary,
    AxisModule,
    AxisResultKind,
    AxisStructureRecurrentSummary,
)


DEFAULT_BATCH_SIZE = 5000
DEFAULT_ITERATOR_CHUNK_SIZE = 5000
DEFAULT_SUMMARY_VERSION = 1


class RecurrentSummaryBuildError(ValueError):
    """Raised when recurrent summary rebuild parameters are invalid."""


def normalize_axis_regulation(value: Any) -> str:
    """
    Normalize an Axis Final regulation label for summary aggregation.

    Empty-like regulation labels are excluded from pattern and dominant-value
    statistics, but their evidence rows still contribute to context_count and
    observation_count.
    """
    value = str(value or "").strip()

    if value.casefold() in {
        "",
        "nan",
        "none",
        "null",
        "<na>",
        "na",
        "n/a",
    }:
        return ""

    return value


def resolve_dominant_regulation(
    regulation_counts: Counter[str],
) -> tuple[str, int]:
    """
    Return a deterministic dominant regulation.

    Ordering:
        1. count descending;
        2. regulation string ascending when counts tie.
    """
    if not regulation_counts:
        return "", 0

    regulation, count = sorted(
        regulation_counts.items(),
        key=lambda item: (-item[1], item[0]),
    )[0]

    return regulation, int(count)


def get_structure_summary_source_queryset() -> QuerySet:
    """
    Aggregate one source row per CanonicalAxis from active context presences.

    AxisContextPresence already represents the union of active artifacts under
    one context, so an Axis present in both Axis Final and Sponge contributes
    only one context occurrence.
    """
    return (
        AxisContextPresence.objects
        .filter(context__is_active=True)
        .values("axis_id")
        .annotate(
            context_count=Count(
                "context_id",
                distinct=True,
            ),
            dataset_count=Count(
                "context__dataset_metadata",
                distinct=True,
            ),
            tcga_dataset_count=Count(
                "context__dataset_metadata",
                filter=Q(
                    context__dataset_source=(
                        AxisDatasetSource.TCGA
                    ),
                ),
                distinct=True,
            ),
            timedb_dataset_count=Count(
                "context__dataset_metadata",
                filter=Q(
                    context__dataset_source=(
                        AxisDatasetSource.TIMEDB
                    ),
                ),
                distinct=True,
            ),
            tcga_context_count=Count(
                "context_id",
                filter=Q(
                    context__dataset_source=(
                        AxisDatasetSource.TCGA
                    ),
                ),
                distinct=True,
            ),
            timedb_context_count=Count(
                "context_id",
                filter=Q(
                    context__dataset_source=(
                        AxisDatasetSource.TIMEDB
                    ),
                ),
                distinct=True,
            ),
            module2_context_count=Count(
                "context_id",
                filter=Q(
                    context__module=AxisModule.MODULE2,
                ),
                distinct=True,
            ),
            module3_context_count=Count(
                "context_id",
                filter=Q(
                    context__module=AxisModule.MODULE3,
                ),
                distinct=True,
            ),
            axis_final_context_count=Count(
                "context_id",
                filter=Q(has_axis_final=True),
                distinct=True,
            ),
            sponge_context_count=Count(
                "context_id",
                filter=Q(has_sponge=True),
                distinct=True,
            ),
            both_result_context_count=Count(
                "context_id",
                filter=Q(
                    has_axis_final=True,
                    has_sponge=True,
                ),
                distinct=True,
            ),
        )
        .order_by("axis_id")
    )


def iter_structure_summary_objects(
    *,
    summary_version: int = DEFAULT_SUMMARY_VERSION,
    iterator_chunk_size: int = DEFAULT_ITERATOR_CHUNK_SIZE,
) -> Iterator[AxisStructureRecurrentSummary]:
    _validate_positive_integer(
        summary_version,
        field_name="summary_version",
    )
    _validate_positive_integer(
        iterator_chunk_size,
        field_name="iterator_chunk_size",
    )

    queryset = get_structure_summary_source_queryset()

    for row in queryset.iterator(
        chunk_size=iterator_chunk_size,
    ):
        yield AxisStructureRecurrentSummary(
            axis_id=row["axis_id"],
            context_count=int(row["context_count"]),
            dataset_count=int(row["dataset_count"]),
            tcga_dataset_count=int(
                row["tcga_dataset_count"]
            ),
            timedb_dataset_count=int(
                row["timedb_dataset_count"]
            ),
            tcga_context_count=int(
                row["tcga_context_count"]
            ),
            timedb_context_count=int(
                row["timedb_context_count"]
            ),
            module2_context_count=int(
                row["module2_context_count"]
            ),
            module3_context_count=int(
                row["module3_context_count"]
            ),
            axis_final_context_count=int(
                row["axis_final_context_count"]
            ),
            sponge_context_count=int(
                row["sponge_context_count"]
            ),
            both_result_context_count=int(
                row["both_result_context_count"]
            ),
            summary_version=summary_version,
        )


def get_structure_summary_source_stats() -> dict:
    source = AxisContextPresence.objects.filter(
        context__is_active=True,
    )

    return {
        "source_presence_count": source.count(),
        "source_axis_count": (
            source.values("axis_id")
            .distinct()
            .count()
        ),
        "source_context_count": (
            source.values("context_id")
            .distinct()
            .count()
        ),
        "source_dataset_count": (
            source.values(
                "context__dataset_metadata"
            )
            .distinct()
            .count()
        ),
    }


@transaction.atomic
def rebuild_axis_structure_recurrent_summary(
    *,
    batch_size: int = DEFAULT_BATCH_SIZE,
    iterator_chunk_size: int | None = None,
    summary_version: int = DEFAULT_SUMMARY_VERSION,
    dry_run: bool = False,
) -> dict:
    """
    Fully rebuild AxisStructureRecurrentSummary from AxisContextPresence.

    Only active contexts are included. Existing summary rows are replaced in
    one transaction so readers do not observe a partially rebuilt table.
    """
    _validate_positive_integer(
        batch_size,
        field_name="batch_size",
    )
    _validate_positive_integer(
        summary_version,
        field_name="summary_version",
    )

    iterator_chunk_size = (
        iterator_chunk_size or batch_size
    )
    _validate_positive_integer(
        iterator_chunk_size,
        field_name="iterator_chunk_size",
    )

    source_stats = get_structure_summary_source_stats()
    current_summary_count = (
        AxisStructureRecurrentSummary.objects.count()
    )

    if dry_run:
        return {
            "success": True,
            "dry_run": True,
            "summary_type": "structure",
            **source_stats,
            "deleted_summary_count": 0,
            "created_summary_count": 0,
            "final_summary_count": current_summary_count,
            "batch_size": batch_size,
            "iterator_chunk_size": iterator_chunk_size,
            "summary_version": summary_version,
        }

    deleted_summary_count, _ = (
        AxisStructureRecurrentSummary.objects
        .all()
        .delete()
    )

    created_summary_count = _bulk_create_stream(
        model=AxisStructureRecurrentSummary,
        objects=iter_structure_summary_objects(
            summary_version=summary_version,
            iterator_chunk_size=iterator_chunk_size,
        ),
        batch_size=batch_size,
    )

    final_summary_count = (
        AxisStructureRecurrentSummary.objects.count()
    )

    if final_summary_count != source_stats["source_axis_count"]:
        raise RuntimeError(
            "AxisStructureRecurrentSummary count does not match "
            "the active presence source axis count: "
            f"summary={final_summary_count}, "
            f"source={source_stats['source_axis_count']}."
        )

    return {
        "success": True,
        "dry_run": False,
        "summary_type": "structure",
        **source_stats,
        "deleted_summary_count": deleted_summary_count,
        "created_summary_count": created_summary_count,
        "final_summary_count": final_summary_count,
        "batch_size": batch_size,
        "iterator_chunk_size": iterator_chunk_size,
        "summary_version": summary_version,
    }


@dataclass(slots=True)
class _AxisFinalAccumulator:
    axis_id: int
    context_ids: set[int] = field(default_factory=set)
    observation_count: int = 0
    regulation_counts: Counter[str] = field(
        default_factory=Counter
    )


def get_axis_final_summary_source_queryset() -> QuerySet:
    """
    Return active Axis Final evidence ordered for streaming aggregation.

    Evidence is used directly rather than AxisContextPresence.has_axis_final so
    regulation statistics cannot become detached from the typed evidence rows.
    """
    return (
        AxisFinalEvidence.objects
        .filter(
            observation__artifact__is_active=True,
            observation__artifact__context__is_active=True,
            observation__artifact__result_kind=(
                AxisResultKind.AXIS_FINAL
            ),
        )
        .order_by(
            "observation__axis_id",
            "observation__artifact__context_id",
            "observation_id",
        )
        .values_list(
            "observation__axis_id",
            "observation__artifact__context_id",
            "axis_regulation",
        )
    )


def _build_axis_final_summary_object(
    accumulator: _AxisFinalAccumulator,
) -> AxisFinalRecurrentSummary:
    dominant_regulation, dominant_count = (
        resolve_dominant_regulation(
            accumulator.regulation_counts
        )
    )
    regulation_pattern_count = len(
        accumulator.regulation_counts
    )

    return AxisFinalRecurrentSummary(
        axis_id=accumulator.axis_id,
        context_count=len(accumulator.context_ids),
        observation_count=(
            accumulator.observation_count
        ),
        regulation_pattern_count=(
            regulation_pattern_count
        ),
        dominant_axis_regulation=(
            dominant_regulation
        ),
        dominant_regulation_count=dominant_count,
        # Consistency is evaluated across non-empty regulation labels.
        # Missing labels do not form an additional regulation pattern.
        regulation_consistent=(
            regulation_pattern_count == 1
        ),
    )


def iter_axis_final_summary_objects(
    *,
    iterator_chunk_size: int = DEFAULT_ITERATOR_CHUNK_SIZE,
) -> Iterator[AxisFinalRecurrentSummary]:
    _validate_positive_integer(
        iterator_chunk_size,
        field_name="iterator_chunk_size",
    )

    current: _AxisFinalAccumulator | None = None

    for (
        axis_id,
        context_id,
        axis_regulation,
    ) in get_axis_final_summary_source_queryset().iterator(
        chunk_size=iterator_chunk_size,
    ):
        axis_id = int(axis_id)
        context_id = int(context_id)

        if current is None or current.axis_id != axis_id:
            if current is not None:
                yield _build_axis_final_summary_object(
                    current
                )

            current = _AxisFinalAccumulator(
                axis_id=axis_id,
            )

        current.context_ids.add(context_id)
        current.observation_count += 1

        regulation = normalize_axis_regulation(
            axis_regulation
        )
        if regulation:
            current.regulation_counts[regulation] += 1

    if current is not None:
        yield _build_axis_final_summary_object(current)


def get_axis_final_summary_source_stats() -> dict:
    source = AxisFinalEvidence.objects.filter(
        observation__artifact__is_active=True,
        observation__artifact__context__is_active=True,
        observation__artifact__result_kind=(
            AxisResultKind.AXIS_FINAL
        ),
    )

    return {
        "source_evidence_count": source.count(),
        "source_axis_count": (
            source.values("observation__axis_id")
            .distinct()
            .count()
        ),
        "source_context_count": (
            source.values(
                "observation__artifact__context_id"
            )
            .distinct()
            .count()
        ),
        "source_nonempty_regulation_count": (
            source.exclude(axis_regulation="").count()
        ),
    }


@transaction.atomic
def rebuild_axis_final_recurrent_summary(
    *,
    batch_size: int = DEFAULT_BATCH_SIZE,
    iterator_chunk_size: int | None = None,
    dry_run: bool = False,
) -> dict:
    """
    Fully rebuild AxisFinalRecurrentSummary from active AxisFinalEvidence.

    context_count counts distinct active contexts containing typed Axis Final
    evidence. observation_count counts all active Axis Final evidence rows.
    Regulation pattern fields ignore empty regulation labels.
    """
    _validate_positive_integer(
        batch_size,
        field_name="batch_size",
    )

    iterator_chunk_size = (
        iterator_chunk_size or batch_size
    )
    _validate_positive_integer(
        iterator_chunk_size,
        field_name="iterator_chunk_size",
    )

    source_stats = get_axis_final_summary_source_stats()
    current_summary_count = (
        AxisFinalRecurrentSummary.objects.count()
    )

    if dry_run:
        return {
            "success": True,
            "dry_run": True,
            "summary_type": "axis_final",
            **source_stats,
            "deleted_summary_count": 0,
            "created_summary_count": 0,
            "final_summary_count": current_summary_count,
            "batch_size": batch_size,
            "iterator_chunk_size": iterator_chunk_size,
        }

    deleted_summary_count, _ = (
        AxisFinalRecurrentSummary.objects
        .all()
        .delete()
    )

    created_summary_count = _bulk_create_stream(
        model=AxisFinalRecurrentSummary,
        objects=iter_axis_final_summary_objects(
            iterator_chunk_size=iterator_chunk_size,
        ),
        batch_size=batch_size,
    )

    final_summary_count = (
        AxisFinalRecurrentSummary.objects.count()
    )

    if final_summary_count != source_stats["source_axis_count"]:
        raise RuntimeError(
            "AxisFinalRecurrentSummary count does not match "
            "the active Axis Final source axis count: "
            f"summary={final_summary_count}, "
            f"source={source_stats['source_axis_count']}."
        )

    return {
        "success": True,
        "dry_run": False,
        "summary_type": "axis_final",
        **source_stats,
        "deleted_summary_count": deleted_summary_count,
        "created_summary_count": created_summary_count,
        "final_summary_count": final_summary_count,
        "batch_size": batch_size,
        "iterator_chunk_size": iterator_chunk_size,
    }


@transaction.atomic
def rebuild_recurrent_axis_summaries(
    *,
    batch_size: int = DEFAULT_BATCH_SIZE,
    iterator_chunk_size: int | None = None,
    summary_version: int = DEFAULT_SUMMARY_VERSION,
    dry_run: bool = False,
) -> dict:
    """
    Rebuild both recurrent summary tables as one logical operation.

    When dry_run=False, the outer transaction guarantees that a failure in the
    Axis Final rebuild also rolls back the structure-summary rebuild.
    """
    structure = rebuild_axis_structure_recurrent_summary(
        batch_size=batch_size,
        iterator_chunk_size=iterator_chunk_size,
        summary_version=summary_version,
        dry_run=dry_run,
    )
    axis_final = rebuild_axis_final_recurrent_summary(
        batch_size=batch_size,
        iterator_chunk_size=iterator_chunk_size,
        dry_run=dry_run,
    )

    consistency = None
    if not dry_run:
        consistency = validate_recurrent_summary_consistency()

    return {
        "success": (
            structure["success"]
            and axis_final["success"]
        ),
        "dry_run": dry_run,
        "structure": structure,
        "axis_final": axis_final,
        "consistency": consistency,
    }


def validate_recurrent_summary_consistency() -> dict:
    """
    Validate the relationship between structure and Axis Final summaries.

    For every Axis with active Axis Final evidence:
        - a structure summary must exist;
        - structure.axis_final_context_count must equal
          axis_final.context_count.

    Any failure indicates stale/corrupt AxisContextPresence or an inconsistent
    active-artifact state. The combined rebuild calls this inside its outer
    transaction, so validation failure rolls back both summary tables.
    """
    axis_final_without_structure = (
        AxisFinalRecurrentSummary.objects
        .filter(axis__recurrent_summary__isnull=True)
        .count()
    )

    structure_without_axis_final = (
        AxisStructureRecurrentSummary.objects
        .filter(
            axis_final_context_count__gt=0,
            axis__axis_final_recurrent_summary__isnull=True,
        )
        .count()
    )

    context_count_mismatch = (
        AxisFinalRecurrentSummary.objects
        .filter(axis__recurrent_summary__isnull=False)
        .exclude(
            context_count=F(
                "axis__recurrent_summary__axis_final_context_count"
            )
        )
        .count()
    )

    if (
        axis_final_without_structure
        or structure_without_axis_final
        or context_count_mismatch
    ):
        raise RuntimeError(
            "Recurrent summary consistency validation failed: "
            f"axis_final_without_structure={axis_final_without_structure}, "
            f"structure_without_axis_final={structure_without_axis_final}, "
            f"context_count_mismatch={context_count_mismatch}. "
            "Rebuild or repair AxisContextPresence before rebuilding "
            "summaries."
        )

    return {
        "axis_final_without_structure": 0,
        "structure_without_axis_final": 0,
        "context_count_mismatch": 0,
    }


def _bulk_create_stream(
    *,
    model,
    objects: Iterable,
    batch_size: int,
) -> int:
    pending = []
    created_count = 0

    for obj in objects:
        pending.append(obj)

        if len(pending) >= batch_size:
            created = model.objects.bulk_create(
                pending,
                batch_size=batch_size,
            )
            created_count += len(created)
            pending = []

    if pending:
        created = model.objects.bulk_create(
            pending,
            batch_size=batch_size,
        )
        created_count += len(created)

    return created_count


def _validate_positive_integer(
    value: int,
    *,
    field_name: str,
) -> None:
    if not isinstance(value, int):
        raise TypeError(
            f"{field_name} must be an integer."
        )

    if value <= 0:
        raise RecurrentSummaryBuildError(
            f"{field_name} must be greater than zero."
        )
