from collections import Counter

from django.db import transaction

from database.models import (
    AxisRecurrentSummary,
    AxisSignatureProjectIndex,
)

from analysis.services.axis_final.axis_signature import (
    parse_axis_signature,
)


DEFAULT_BATCH_SIZE = 5000


class AxisRecurrentSummaryBuildError(ValueError):
    pass


def normalize_regulation(value) -> str:
    value = str(value or "").strip()

    if value.lower() in {
        "",
        "nan",
        "none",
        "null",
        "<na>",
    }:
        return ""

    return value


def create_empty_accumulator(
    *,
    axis_signature: str,
) -> dict:
    axis_parts = parse_axis_signature(axis_signature)

    return {
        "axis_signature": axis_signature,

        "axis_type": axis_parts["axis_type"],
        "miRNA": axis_parts["miRNA"],
        "mRNA": axis_parts["mRNA"],
        "lncRNA": axis_parts["lncRNA"],
        "circRNA": axis_parts["circRNA"],

        "project_ids": set(),
        "dataset_names": set(),

        "tcga_project_ids": set(),
        "timedb_project_ids": set(),

        "regulation_counts": Counter(),
    }


def add_index_row_to_accumulator(
    *,
    accumulator: dict,
    row: dict,
) -> None:
    project_id = row["project_id"]
    source = str(row.get("source") or "").strip().upper()
    dataset_name = str(row.get("dataset_name") or "").strip()
    axis_regulation = normalize_regulation(
        row.get("axis_regulation")
    )

    accumulator["project_ids"].add(project_id)

    if dataset_name:
        accumulator["dataset_names"].add(dataset_name)

    if source == "TCGA":
        accumulator["tcga_project_ids"].add(project_id)

    elif source == "TIMEDB":
        accumulator["timedb_project_ids"].add(project_id)

    if axis_regulation:
        accumulator["regulation_counts"][axis_regulation] += 1


def resolve_dominant_regulation(
    regulation_counts: Counter,
) -> tuple[str, int]:
    if not regulation_counts:
        return "", 0

    # First sort by:
    # 1. count descending
    # 2. regulation name ascending
    #
    # This ensures deterministic results when counts tie.
    dominant_regulation, dominant_count = sorted(
        regulation_counts.items(),
        key=lambda item: (
            -item[1],
            item[0],
        ),
    )[0]

    return dominant_regulation, int(dominant_count)


def build_summary_object_from_accumulator(
    accumulator: dict,
) -> AxisRecurrentSummary:
    regulation_counts = accumulator["regulation_counts"]

    dominant_regulation, dominant_count = (
        resolve_dominant_regulation(regulation_counts)
    )

    regulation_pattern_count = len(regulation_counts)

    return AxisRecurrentSummary(
        axis_signature=accumulator["axis_signature"],

        axis_type=accumulator["axis_type"],
        miRNA=accumulator["miRNA"],
        mRNA=accumulator["mRNA"],
        lncRNA=accumulator["lncRNA"],
        circRNA=accumulator["circRNA"],

        project_count=len(accumulator["project_ids"]),
        dataset_count=len(accumulator["dataset_names"]),

        tcga_project_count=len(
            accumulator["tcga_project_ids"]
        ),
        timedb_project_count=len(
            accumulator["timedb_project_ids"]
        ),

        regulation_pattern_count=regulation_pattern_count,
        dominant_axis_regulation=dominant_regulation,
        dominant_regulation_count=dominant_count,

        # Empty regulation data does not count as consistent.
        regulation_consistent=(
            regulation_pattern_count == 1
        ),
    )


def get_recurrent_summary_source_queryset(
    *,
    active_projects_only: bool = True,
):
    """
    AxisSignatureProjectIndex has one row per:
        axis_signature + project

    Therefore each row represents one project occurrence and can be used
    directly for project-level recurrent statistics.
    """
    qs = (
        AxisSignatureProjectIndex.objects
        .exclude(axis_signature="")
        .order_by(
            "axis_signature",
            "project_id",
        )
        .values(
            "axis_signature",
            "project_id",
            "source",
            "dataset_name",
            "axis_regulation",
        )
    )

    if active_projects_only:
        qs = qs.filter(project__is_active=True)

    return qs


def iter_axis_recurrent_summary_objects(
    *,
    active_projects_only: bool = True,
    iterator_chunk_size: int = DEFAULT_BATCH_SIZE,
):
    """
    Stream AxisSignatureProjectIndex ordered by axis_signature and emit
    one AxisRecurrentSummary object per signature.

    This avoids loading the whole index table into memory.
    """
    if iterator_chunk_size <= 0:
        raise AxisRecurrentSummaryBuildError(
            "iterator_chunk_size must be greater than 0."
        )

    qs = get_recurrent_summary_source_queryset(
        active_projects_only=active_projects_only,
    )

    current_signature = None
    accumulator = None

    for row in qs.iterator(
        chunk_size=iterator_chunk_size,
    ):
        axis_signature = str(
            row.get("axis_signature") or ""
        ).strip()

        if not axis_signature:
            continue

        if current_signature != axis_signature:
            if accumulator is not None:
                yield build_summary_object_from_accumulator(
                    accumulator
                )

            current_signature = axis_signature
            accumulator = create_empty_accumulator(
                axis_signature=axis_signature,
            )

        add_index_row_to_accumulator(
            accumulator=accumulator,
            row=row,
        )

    if accumulator is not None:
        yield build_summary_object_from_accumulator(
            accumulator
        )


def get_axis_recurrent_summary_source_stats(
    *,
    active_projects_only: bool = True,
) -> dict:
    qs = get_recurrent_summary_source_queryset(
        active_projects_only=active_projects_only,
    )

    return {
        "source_index_row_count": qs.count(),
        "source_axis_signature_count": (
            qs.values("axis_signature")
            .distinct()
            .count()
        ),
    }


@transaction.atomic
def rebuild_axis_recurrent_summary(
    *,
    batch_size: int = DEFAULT_BATCH_SIZE,
    iterator_chunk_size: int | None = None,
    active_projects_only: bool = True,
    dry_run: bool = False,
) -> dict:
    """
    Fully rebuild AxisRecurrentSummary from AxisSignatureProjectIndex.

    Recommended order:
        1. import Module2
        2. import Module3
        3. rebuild AxisSignatureProjectIndex
        4. rebuild AxisRecurrentSummary
    """
    if batch_size <= 0:
        raise AxisRecurrentSummaryBuildError(
            "batch_size must be greater than 0."
        )

    iterator_chunk_size = (
        iterator_chunk_size or batch_size
    )

    if iterator_chunk_size <= 0:
        raise AxisRecurrentSummaryBuildError(
            "iterator_chunk_size must be greater than 0."
        )

    source_stats = get_axis_recurrent_summary_source_stats(
        active_projects_only=active_projects_only,
    )

    if dry_run:
        return {
            "success": True,
            "dry_run": True,

            **source_stats,

            "deleted_summary_count": 0,
            "created_summary_count": 0,
            "final_summary_count": (
                AxisRecurrentSummary.objects.count()
            ),

            "batch_size": batch_size,
            "iterator_chunk_size": iterator_chunk_size,
            "active_projects_only": active_projects_only,
        }

    deleted_summary_count, _ = (
        AxisRecurrentSummary.objects.all().delete()
    )

    pending_objects = []
    created_summary_count = 0

    for summary_object in iter_axis_recurrent_summary_objects(
        active_projects_only=active_projects_only,
        iterator_chunk_size=iterator_chunk_size,
    ):
        pending_objects.append(summary_object)

        if len(pending_objects) >= batch_size:
            created_objects = (
                AxisRecurrentSummary.objects.bulk_create(
                    pending_objects,
                    batch_size=batch_size,
                )
            )

            created_summary_count += len(created_objects)
            pending_objects = []

    if pending_objects:
        created_objects = (
            AxisRecurrentSummary.objects.bulk_create(
                pending_objects,
                batch_size=batch_size,
            )
        )

        created_summary_count += len(created_objects)

    final_summary_count = (
        AxisRecurrentSummary.objects.count()
    )

    return {
        "success": True,
        "dry_run": False,

        **source_stats,

        "deleted_summary_count": deleted_summary_count,
        "created_summary_count": created_summary_count,
        "final_summary_count": final_summary_count,

        "batch_size": batch_size,
        "iterator_chunk_size": iterator_chunk_size,
        "active_projects_only": active_projects_only,
    }
