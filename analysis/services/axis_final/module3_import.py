import hashlib
import math
from pathlib import Path

import pandas as pd
from django.db import transaction

from analysis.utils.workflow_detail_utils.workflow_axis_final_utils import (
    HYBRID_REFERENCE_AXIS_FINAL_COLUMNS,
    HYBRID_REFERENCE_AXIS_FINAL_REQUIRED_COLUMNS,
    WORKFLOW_AXIS_FINAL_FILENAME_SUFFIX,
    read_axis_final_file_by_path,
    normalize_workflow_axis_final_dataframe,
)

from database.models import (
    DatasetAxisFinalProject,
    DatasetAxisFinalOccurrence,
)

from database.utils.dataset_annotation_utils.path_utils import (
    get_timedb_json_group_by_fields,
)

from analysis.services.axis_final.axis_signature import (
    add_axis_signature_to_dataframe,
)


MODULE3_SOURCE = DatasetAxisFinalProject.Source.TIMEDB
MODULE3_MODULE = DatasetAxisFinalProject.Module.MODULE3

MODULE3_COMMON_GROUP_TYPE = DatasetAxisFinalProject.GroupType.COMMON
MODULE3_GRADE_GROUP_TYPE = DatasetAxisFinalProject.GroupType.GRADE
MODULE3_STAGE_GROUP_TYPE = DatasetAxisFinalProject.GroupType.STAGE


class Module3AxisFinalImportError(ValueError):
    pass


def safe_str(value) -> str:
    if value is None:
        return ""

    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass

    return str(value).strip()


def safe_float(value):
    if value is None:
        return None

    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass

    try:
        value = float(value)
    except (TypeError, ValueError):
        return None

    if not math.isfinite(value):
        return None

    return value


def calculate_file_sha256(file_path: Path) -> str:
    hasher = hashlib.sha256()

    with Path(file_path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            hasher.update(chunk)

    return hasher.hexdigest()


def iter_module3_annotation_dirs(module3_root_dir) -> list[Path]:
    module3_root_dir = Path(module3_root_dir).resolve()

    if not module3_root_dir.exists() or not module3_root_dir.is_dir():
        raise Module3AxisFinalImportError(
            f"Module3 result directory is not available: {module3_root_dir}"
        )

    return sorted(
        item
        for item in module3_root_dir.iterdir()
        if item.is_dir() and not item.name.startswith(".")
    )


def parse_module3_dir_name(annotation_dir: Path) -> dict:
    """
    Parse Module3 annotation directory.

    Examples:
        GSE20194
            dataset_name = GSE20194
            group_type = common

        GSE20194_grade
            dataset_name = GSE20194
            group_type = grade

        GSE20194_stage
            dataset_name = GSE20194
            group_type = stage
    """
    dir_name = annotation_dir.name

    if dir_name.endswith("_grade"):
        dataset_name = dir_name[: -len("_grade")]

        return {
            "annotation_dir_name": dir_name,
            "dataset_name": dataset_name,
            "annotation_file_prefix": dataset_name,
            "group_type": MODULE3_GRADE_GROUP_TYPE,
            "group_by_values": ["Grade"],
        }

    if dir_name.endswith("_stage"):
        dataset_name = dir_name[: -len("_stage")]

        return {
            "annotation_dir_name": dir_name,
            "dataset_name": dataset_name,
            "annotation_file_prefix": dataset_name,
            "group_type": MODULE3_STAGE_GROUP_TYPE,
            "group_by_values": ["Stage"],
        }

    dataset_name = dir_name

    group_by_values = [
        str(item).strip()
        for item in get_timedb_json_group_by_fields(dataset_name)
        if str(item or "").strip()
    ]

    if not group_by_values:
        group_by_values = [""]

    return {
        "annotation_dir_name": dir_name,
        "dataset_name": dataset_name,
        "annotation_file_prefix": dataset_name,
        "group_type": MODULE3_COMMON_GROUP_TYPE,
        "group_by_values": group_by_values,
    }


def build_module3_axis_final_file(
    *,
    annotation_dir: Path,
    dataset_name: str,
) -> Path:
    """
    File naming rule:
        {dataset_name}_ceRNA_axis_final.csv

    Examples:
        GSE20194/GSE20194_ceRNA_axis_final.csv
        GSE20194_grade/GSE20194_ceRNA_axis_final.csv
        GSE20194_stage/GSE20194_ceRNA_axis_final.csv
    """
    return (
        annotation_dir / f"{dataset_name}{WORKFLOW_AXIS_FINAL_FILENAME_SUFFIX}"
    ).resolve()


def read_module3_axis_final_dataframe(axis_final_file: Path) -> pd.DataFrame:
    _, df = read_axis_final_file_by_path(
        file_path=axis_final_file,
        required_columns=HYBRID_REFERENCE_AXIS_FINAL_REQUIRED_COLUMNS,
    )

    df = normalize_workflow_axis_final_dataframe(
        df=df,
        columns=HYBRID_REFERENCE_AXIS_FINAL_COLUMNS,
    )

    df = add_axis_signature_to_dataframe(df)

    return df


@transaction.atomic
def import_one_module3_axis_final_project(
    *,
    annotation_dir: Path,
    group_by: str,
    dry_run: bool = False,
) -> dict:
    context = parse_module3_dir_name(annotation_dir)

    dataset_name = context["dataset_name"]
    group_type = context["group_type"]
    annotation_dir_name = context["annotation_dir_name"]
    annotation_file_prefix = context["annotation_file_prefix"]

    group_by = safe_str(group_by)

    axis_final_file = build_module3_axis_final_file(
        annotation_dir=annotation_dir,
        dataset_name=dataset_name,
    )

    if not axis_final_file.exists() or not axis_final_file.is_file():
        return {
            "success": False,
            "skipped": True,
            "reason": "axis_final_file_not_found",
            "annotation_dir_name": annotation_dir_name,
            "dataset_name": dataset_name,
            "group_type": group_type,
            "group_by": group_by,
            "axis_final_file": str(axis_final_file),
        }

    df = read_module3_axis_final_dataframe(axis_final_file)

    file_sha256 = calculate_file_sha256(axis_final_file)
    row_count = int(df.shape[0])

    if dry_run:
        return {
            "success": True,
            "dry_run": True,
            "annotation_dir_name": annotation_dir_name,
            "dataset_name": dataset_name,
            "group_type": group_type,
            "group_by": group_by,
            "axis_final_file": str(axis_final_file),
            "file_sha256": file_sha256,
            "row_count": row_count,
        }

    project, _ = DatasetAxisFinalProject.objects.update_or_create(
        source=MODULE3_SOURCE,
        module=MODULE3_MODULE,
        dataset_name=dataset_name,
        group_type=group_type,
        group_by=group_by,
        defaults={
            "annotation_dir_name": annotation_dir_name,
            "annotation_file_prefix": annotation_file_prefix,
            "axis_final_file_name": axis_final_file.name,
            "axis_final_file_path": str(axis_final_file),
            "file_sha256": file_sha256,
            "row_count": row_count,
            "is_active": True,
        },
    )

    project.axis_occurrences.all().delete()

    occurrences = []

    for row_index, row in df.iterrows():
        occurrences.append(
            DatasetAxisFinalOccurrence(
                project=project,
                row_index=int(row_index),

                axis_signature=safe_str(row.get("axis_signature")),
                axis_id=safe_str(row.get("axis_id")),
                axis_type=safe_str(row.get("axis_type")),
                axis_regulation=safe_str(row.get("axis_regulation")),

                miRNA=safe_str(row.get("miRNA")),
                mRNA=safe_str(row.get("mRNA")),
                lncRNA=safe_str(row.get("lncRNA")),
                circRNA=safe_str(row.get("circRNA")),

                mRNA_log2FC=safe_float(row.get("mRNA_log2FC")),
                mRNA_regulation=safe_str(row.get("mRNA_regulation")),

                miRNA_log2FC=safe_float(row.get("miRNA_log2FC")),
                miRNA_regulation=safe_str(row.get("miRNA_regulation")),

                lncRNA_log2FC=safe_float(row.get("lncRNA_log2FC")),
                lncRNA_regulation=safe_str(row.get("lncRNA_regulation")),

                circRNA_log2FC=safe_float(row.get("circRNA_log2FC")),
                circRNA_regulation=safe_str(row.get("circRNA_regulation")),
            )
        )

    DatasetAxisFinalOccurrence.objects.bulk_create(
        occurrences,
        batch_size=1000,
    )

    return {
        "success": True,
        "skipped": False,
        "annotation_dir_name": annotation_dir_name,
        "dataset_name": dataset_name,
        "group_type": group_type,
        "group_by": group_by,
        "project_id": project.id,
        "axis_final_file": str(axis_final_file),
        "file_sha256": file_sha256,
        "row_count": row_count,
        "imported_occurrence_count": len(occurrences),
    }


def import_one_module3_annotation_dir(
    *,
    annotation_dir: Path,
    dry_run: bool = False,
) -> list[dict]:
    """
    Import one Module3 directory.

    For common:
        if JSON contains multiple group_by fields, this function creates
        one DatasetAxisFinalProject per group_by value.

    For grade / stage:
        creates exactly one DatasetAxisFinalProject.
    """
    context = parse_module3_dir_name(annotation_dir)

    results = []

    for group_by in context["group_by_values"]:
        try:
            result = import_one_module3_axis_final_project(
                annotation_dir=annotation_dir,
                group_by=group_by,
                dry_run=dry_run,
            )
        except Exception as e:
            result = {
                "success": False,
                "skipped": False,
                "annotation_dir_name": annotation_dir.name,
                "dataset_name": context.get("dataset_name"),
                "group_type": context.get("group_type"),
                "group_by": group_by,
                "error": str(e),
            }

        results.append(result)

    return results


def import_module3_axis_final_projects(
    *,
    module3_root_dir,
    dry_run: bool = False,
) -> dict:
    annotation_dirs = iter_module3_annotation_dirs(module3_root_dir)

    results = []

    for annotation_dir in annotation_dirs:
        dir_results = import_one_module3_annotation_dir(
            annotation_dir=annotation_dir,
            dry_run=dry_run,
        )

        results.extend(dir_results)

    imported_count = sum(
        1
        for item in results
        if (
            item.get("success")
            and not item.get("skipped")
            and not item.get("dry_run")
        )
    )

    dry_run_count = sum(
        1
        for item in results
        if item.get("success") and item.get("dry_run")
    )

    skipped_count = sum(
        1
        for item in results
        if item.get("skipped")
    )

    failed_count = sum(
        1
        for item in results
        if not item.get("success") and not item.get("skipped")
    )

    return {
        "success": failed_count == 0,
        "source": MODULE3_SOURCE,
        "module": MODULE3_MODULE,
        "module3_root_dir": str(Path(module3_root_dir).resolve()),
        "annotation_dir_count": len(annotation_dirs),
        "project_count": len(results),
        "imported_count": imported_count,
        "dry_run_count": dry_run_count,
        "skipped_count": skipped_count,
        "failed_count": failed_count,
        "results": results,
    }
