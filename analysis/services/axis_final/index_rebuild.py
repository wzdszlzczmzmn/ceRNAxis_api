from django.db import transaction

from database.models import (
    DatasetAxisFinalOccurrence,
    AxisSignatureProjectIndex,
)


DEFAULT_BATCH_SIZE = 5000


class AxisSignatureProjectIndexRebuildError(ValueError):
    pass


def build_index_row_from_occurrence(occurrence):
    """
    Build one AxisSignatureProjectIndex object from one
    DatasetAxisFinalOccurrence object.

    Source of truth:
        DatasetAxisFinalOccurrence
        DatasetAxisFinalProject

    Rebuild target:
        AxisSignatureProjectIndex
    """
    project = occurrence.project

    return AxisSignatureProjectIndex(
        axis_signature=occurrence.axis_signature,

        project=project,
        occurrence=occurrence,

        source=project.source,
        module=project.module,
        dataset_name=project.dataset_name,
        group_type=project.group_type,
        group_by=project.group_by,

        axis_id=occurrence.axis_id,
        axis_type=occurrence.axis_type,
        axis_regulation=occurrence.axis_regulation,
    )


def get_axis_signature_index_source_queryset(
    *,
    active_projects_only: bool = True,
):
    """
    Return occurrence queryset used to rebuild AxisSignatureProjectIndex.
    """
    qs = (
        DatasetAxisFinalOccurrence.objects
        .select_related("project")
        .exclude(axis_signature="")
        .order_by("id")
    )

    if active_projects_only:
        qs = qs.filter(project__is_active=True)

    return qs


@transaction.atomic
def clear_axis_signature_project_index():
    deleted_count, _ = AxisSignatureProjectIndex.objects.all().delete()

    return deleted_count


def rebuild_axis_signature_project_index(
    *,
    batch_size: int = DEFAULT_BATCH_SIZE,
    dry_run: bool = False,
    active_projects_only: bool = True,
    clear_existing: bool = True,
) -> dict:
    """
    Rebuild AxisSignatureProjectIndex from DatasetAxisFinalOccurrence.

    Recommended usage:
        after importing all Module2 / Module3 axis_final reference data.

    Args:
        batch_size:
            bulk_create batch size.
        dry_run:
            If True, only count rows and validate source queryset.
        active_projects_only:
            If True, only index occurrences whose project is active.
        clear_existing:
            If True, clear old index rows before rebuilding.

    Returns:
        summary dict.
    """
    if batch_size <= 0:
        raise AxisSignatureProjectIndexRebuildError(
            "batch_size must be greater than 0."
        )

    source_qs = get_axis_signature_index_source_queryset(
        active_projects_only=active_projects_only,
    )

    source_count = source_qs.count()

    if dry_run:
        return {
            "success": True,
            "dry_run": True,
            "source_occurrence_count": source_count,
            "deleted_index_count": 0,
            "created_index_count": 0,
            "batch_size": batch_size,
            "active_projects_only": active_projects_only,
            "clear_existing": clear_existing,
        }

    deleted_index_count = 0

    if clear_existing:
        deleted_index_count = clear_axis_signature_project_index()

    rows = []
    created_index_count = 0

    for occurrence in source_qs.iterator(chunk_size=batch_size):
        rows.append(
            build_index_row_from_occurrence(occurrence)
        )

        if len(rows) >= batch_size:
            created = AxisSignatureProjectIndex.objects.bulk_create(
                rows,
                batch_size=batch_size,
                ignore_conflicts=True,
            )

            created_index_count += len(created)
            rows = []

    if rows:
        created = AxisSignatureProjectIndex.objects.bulk_create(
            rows,
            batch_size=batch_size,
            ignore_conflicts=True,
        )

        created_index_count += len(created)

    final_index_count = AxisSignatureProjectIndex.objects.count()

    return {
        "success": True,
        "dry_run": False,
        "source_occurrence_count": source_count,
        "deleted_index_count": deleted_index_count,
        "created_index_count": created_index_count,
        "final_index_count": final_index_count,
        "batch_size": batch_size,
        "active_projects_only": active_projects_only,
        "clear_existing": clear_existing,
    }
