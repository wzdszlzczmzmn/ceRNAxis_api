from dataclasses import dataclass
from pathlib import Path
import re

from analysis.utils.custom_list_query_task_utils import (
    validate_task_name_for_filename as validate_custom_task_name, get_immune_result_file_name,
    get_immune_result_file_path,
)
from analysis.utils.paired_cohort_task_utils import (
    get_paired_cohort_task_output_dir,
    validate_task_name_for_filename as validate_paired_cohort_task_name,
)


@dataclass(frozen=True)
class DownloadableResultFile:
    path: Path
    arcname: str
    required: bool = True


PAIRED_COHORT_DEG_RNA_TYPES = [
    "mRNA",
    "miRNA",
    "lncRNA",
]


def sanitize_archive_name_part(value: str) -> str:
    value = str(value or "").strip()

    if not value:
        return "task"

    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = value.strip("._-")

    return value or "task"


def get_short_uuid(task) -> str:
    return str(task.uuid).split("-")[0]


def build_archive_root_dir(task) -> str:
    safe_task_name = sanitize_archive_name_part(task.task_name)
    short_uuid = get_short_uuid(task)

    return f"{safe_task_name}_{short_uuid}_result"


def build_custom_list_query_archive_name(task) -> str:
    return f"{build_archive_root_dir(task)}.zip"


def build_paired_cohort_archive_name(task) -> str:
    return f"{build_archive_root_dir(task)}.zip"


def resolve_custom_list_query_result_files(task) -> list[DownloadableResultFile]:
    """
    CustomListQueryTask result files:

        output/{task_name}_map_immune_gene_axis.csv
    """

    validate_custom_task_name(task.task_name)

    archive_root = build_archive_root_dir(task)

    result_filename = get_immune_result_file_name(task)
    result_file_path = get_immune_result_file_path(task)

    return [
        DownloadableResultFile(
            path=result_file_path,
            arcname=f"{archive_root}/{result_filename}",
            required=True,
        )
    ]


def resolve_paired_cohort_result_files(task) -> list[DownloadableResultFile]:
    """
    PairedCohortTask result files:

        output/{task_name}_ceRNA_axis.csv
        output/{task_name}_ceRNA_background.csv
        output/{task_name}_ceRNA_corr.csv
        output/{task_name}_map_immune_axis.csv
        output/{task_name}_ceRNA_network.csv
        output/{task_name}_{deg_method}_{rna_type}.csv

    Notes:
        {task_name}_ceRNA_network.csv is treated as optional here.
    """

    validate_paired_cohort_task_name(task.task_name)

    task_name = str(task.task_name).strip()
    deg_method = str(task.deg_method).strip()

    if deg_method not in ["limma", "deseq2"]:
        raise ValueError(
            "Invalid deg_method. Allowed values are: limma, deseq2."
        )

    output_dir = get_paired_cohort_task_output_dir(task)
    archive_root = build_archive_root_dir(task)

    result_files = [
        DownloadableResultFile(
            path=(output_dir / f"{task_name}_ceRNA_axis.csv").resolve(),
            arcname=f"{archive_root}/{task_name}_ceRNA_axis.csv",
            required=True,
        ),
        DownloadableResultFile(
            path=(output_dir / f"{task_name}_ceRNA_background.csv").resolve(),
            arcname=f"{archive_root}/{task_name}_ceRNA_background.csv",
            required=True,
        ),
        DownloadableResultFile(
            path=(output_dir / f"{task_name}_ceRNA_corr.csv").resolve(),
            arcname=f"{archive_root}/{task_name}_ceRNA_corr.csv",
            required=True,
        ),
        DownloadableResultFile(
            path=(output_dir / f"{task_name}_map_immune_axis.csv").resolve(),
            arcname=f"{archive_root}/{task_name}_map_immune_axis.csv",
            required=True,
        ),
        DownloadableResultFile(
            path=(output_dir / f"{task_name}_ceRNA_network.csv").resolve(),
            arcname=f"{archive_root}/{task_name}_ceRNA_network.csv",
            required=False,
        ),
    ]

    for rna_type in PAIRED_COHORT_DEG_RNA_TYPES:
        filename = f"{task_name}_{deg_method}_{rna_type}.csv"

        result_files.append(
            DownloadableResultFile(
                path=(output_dir / filename).resolve(),
                arcname=f"{archive_root}/{filename}",
                required=True,
            )
        )

    return result_files
