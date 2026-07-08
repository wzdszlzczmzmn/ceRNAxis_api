import hashlib
import math
from pathlib import Path

import pandas as pd
from django.db import transaction

from analysis.utils.workflow_detail_utils.workflow_axis_final_utils import (
    PAIRED_COHORT_AXIS_FINAL_COLUMNS,
    PAIRED_COHORT_AXIS_FINAL_REQUIRED_COLUMNS,
    WORKFLOW_AXIS_FINAL_FILENAME_SUFFIX,
    read_axis_final_file_by_path,
    normalize_workflow_axis_final_dataframe,
)

from database.models import (
    DatasetAxisFinalProject,
    DatasetAxisFinalOccurrence,
)

from analysis.services.axis_final.axis_signature import (
    add_axis_signature_to_dataframe,
)


MODULE2_SOURCE = DatasetAxisFinalProject.Source.TCGA
MODULE2_MODULE = DatasetAxisFinalProject.Module.MODULE2
MODULE2_GROUP_TYPE = DatasetAxisFinalProject.GroupType.NONE
MODULE2_GROUP_BY = ""


class Module2AxisFinalImportError(ValueError):
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


def iter_module2_project_dirs(module2_root_dir) -> list[Path]:
    module2_root_dir = Path(module2_root_dir).resolve()

    if not module2_root_dir.exists() or not module2_root_dir.is_dir():
        raise Module2AxisFinalImportError(
            f"Module2 result directory is not available: {module2_root_dir}"
        )

    return sorted(
        item
        for item in module2_root_dir.iterdir()
        if item.is_dir() and not item.name.startswith(".")
    )


def build_module2_project_context(project_dir: Path) -> dict:
    """
    Example:
        project_dir.name = TCGA_ACC

    Then:
        dataset_name = TCGA_ACC_mRNA
        annotation_dir_name = TCGA_ACC
        annotation_file_prefix = TCGA_ACC
        axis_final_file = TCGA_ACC/TCGA_ACC_ceRNA_axis_final.csv
    """
    folder_name = project_dir.name

    dataset_name = f"{folder_name}_mRNA"
    annotation_dir_name = folder_name
    annotation_file_prefix = folder_name

    axis_final_file = (
        project_dir / f"{folder_name}{WORKFLOW_AXIS_FINAL_FILENAME_SUFFIX}"
    ).resolve()

    return {
        "folder_name": folder_name,
        "dataset_name": dataset_name,
        "annotation_dir_name": annotation_dir_name,
        "annotation_file_prefix": annotation_file_prefix,
        "axis_final_file": axis_final_file,
    }


def read_module2_axis_final_dataframe(axis_final_file: Path) -> pd.DataFrame:
    _, df = read_axis_final_file_by_path(
        file_path=axis_final_file,
        required_columns=PAIRED_COHORT_AXIS_FINAL_REQUIRED_COLUMNS,
    )

    df = normalize_workflow_axis_final_dataframe(
        df=df,
        columns=PAIRED_COHORT_AXIS_FINAL_COLUMNS,
    )

    df = add_axis_signature_to_dataframe(df)

    return df


@transaction.atomic
def import_one_module2_axis_final_project(
    *,
    project_dir: Path,
    dry_run: bool = False,
) -> dict:
    context = build_module2_project_context(project_dir)

    folder_name = context["folder_name"]
    dataset_name = context["dataset_name"]
    annotation_dir_name = context["annotation_dir_name"]
    annotation_file_prefix = context["annotation_file_prefix"]
    axis_final_file = context["axis_final_file"]

    if not axis_final_file.exists() or not axis_final_file.is_file():
        return {
            "success": False,
            "skipped": True,
            "reason": "axis_final_file_not_found",
            "folder_name": folder_name,
            "dataset_name": dataset_name,
            "axis_final_file": str(axis_final_file),
        }

    df = read_module2_axis_final_dataframe(axis_final_file)

    file_sha256 = calculate_file_sha256(axis_final_file)
    row_count = int(df.shape[0])

    if dry_run:
        return {
            "success": True,
            "dry_run": True,
            "folder_name": folder_name,
            "dataset_name": dataset_name,
            "axis_final_file": str(axis_final_file),
            "file_sha256": file_sha256,
            "row_count": row_count,
        }

    project, _ = DatasetAxisFinalProject.objects.update_or_create(
        source=MODULE2_SOURCE,
        module=MODULE2_MODULE,
        dataset_name=dataset_name,
        group_type=MODULE2_GROUP_TYPE,
        group_by=MODULE2_GROUP_BY,
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

    # 重新导入时覆盖该 dataset_name 对应的旧 axis occurrences。
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
        "folder_name": folder_name,
        "dataset_name": dataset_name,
        "project_id": project.id,
        "axis_final_file": str(axis_final_file),
        "file_sha256": file_sha256,
        "row_count": row_count,
        "imported_occurrence_count": len(occurrences),
    }


def import_module2_axis_final_projects(
    *,
    module2_root_dir,
    dry_run: bool = False,
) -> dict:
    project_dirs = iter_module2_project_dirs(module2_root_dir)

    results = []

    for project_dir in project_dirs:
        try:
            result = import_one_module2_axis_final_project(
                project_dir=project_dir,
                dry_run=dry_run,
            )
        except Exception as e:
            result = {
                "success": False,
                "skipped": False,
                "folder_name": project_dir.name,
                "error": str(e),
            }

        results.append(result)

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
        "source": MODULE2_SOURCE,
        "module": MODULE2_MODULE,
        "group_type": MODULE2_GROUP_TYPE,
        "group_by": MODULE2_GROUP_BY,
        "module2_root_dir": str(Path(module2_root_dir).resolve()),
        "project_dir_count": len(project_dirs),
        "imported_count": imported_count,
        "dry_run_count": dry_run_count,
        "skipped_count": skipped_count,
        "failed_count": failed_count,
        "results": results,
    }
