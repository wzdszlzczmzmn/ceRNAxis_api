from dataclasses import dataclass
from pathlib import Path
import re

from analysis.utils.custom_list_query_task_utils import (
    validate_task_name_for_filename as validate_custom_task_name, get_immune_result_file_name,
    get_immune_result_file_path,
)
from analysis.utils.hybrid_reference_task_utils import get_hybrid_reference_task_output_dir
from analysis.utils.paired_cohort_task_utils import (
    get_paired_cohort_task_output_dir,
    validate_task_name_for_filename as validate_paired_cohort_task_name,
)


@dataclass(frozen=True)
class DownloadableResultFile:
    path: Path
    arcname: str
    required: bool = True
    is_directory: bool = False


PAIRED_COHORT_DEG_RNA_TYPES = [
    "mRNA",
    "miRNA",
    "lncRNA",
    "circRNA",
]


HYBRID_REFERENCE_DEG_RNA_TYPES = [
    "mRNA",
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


def build_hybrid_reference_archive_name(task) -> str:
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
    PairedCohortTask result files.

    Download rule:
        - All configured result files are optional.
        - Existing files are packed into the zip archive.
        - Missing files are skipped.
        - {task_name}_CMdrug_result/ must be included as a directory entry,
          even if the directory is empty.
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

    paired_cohort_result_filenames = [
        f"{task_name}_ceRNA_axis.csv",
        f"{task_name}_ceRNA_axis_final.csv",
        f"{task_name}_ceRNA_background.csv",
        f"{task_name}_ceRNA_corr.csv",
        f"{task_name}_ceRNA_network.csv",
        f"{task_name}_CMap.csv",
        f"{task_name}_map_immune_axis.csv",
        f"{task_name}_mRNA_gsea.csv",
        f"{task_name}_survival_analysis.csv",
    ]

    result_files = [
        DownloadableResultFile(
            path=(output_dir / filename).resolve(),
            arcname=f"{archive_root}/{filename}",
            required=False,
        )
        for filename in paired_cohort_result_filenames
    ]

    for rna_type in PAIRED_COHORT_DEG_RNA_TYPES:
        filename = f"{task_name}_{deg_method}_{rna_type}.csv"

        result_files.append(
            DownloadableResultFile(
                path=(output_dir / filename).resolve(),
                arcname=f"{archive_root}/{filename}",
                required=False,
            )
        )

    cm_drug_result_dirname = f"{task_name}_CMdrug_result"

    result_files.append(
        DownloadableResultFile(
            path=(output_dir / cm_drug_result_dirname).resolve(),
            arcname=f"{archive_root}/{cm_drug_result_dirname}/",
            required=False,
            is_directory=True,
        )
    )

    return result_files


def resolve_hybrid_reference_result_files(task) -> list[DownloadableResultFile]:
    """
    HybridReferenceTask result files.

    Download rule:
        - All configured result files are optional.
        - Existing files are packed into the zip archive.
        - Missing files are skipped.
        - {task_name}_CMdrug_result/ must be included as a directory entry,
          even if the physical directory does not exist or is empty.

    Hybrid Reference DEG files:
        - {task_name}_{deg_method}_mRNA.csv
        - {task_name}_{deg_method}_mRNA_intersect.csv
        - {task_name}_{deg_method}_mRNA_venn.csv
    """

    validate_paired_cohort_task_name(task.task_name)

    task_name = str(task.task_name).strip()
    deg_method = str(task.deg_method).strip()

    if deg_method not in ["limma", "deseq2"]:
        raise ValueError(
            "Invalid deg_method. Allowed values are: limma, deseq2."
        )

    output_dir = get_hybrid_reference_task_output_dir(task)
    archive_root = build_archive_root_dir(task)

    hybrid_reference_result_filenames = [
        f"{task_name}_ceRNA_axis.csv",
        f"{task_name}_ceRNA_axis_final.csv",
        f"{task_name}_ceRNA_background.csv",
        f"{task_name}_ceRNA_corr.csv",
        f"{task_name}_ceRNA_network.csv",
        f"{task_name}_CMap.csv",
        f"{task_name}_map_immune_axis.csv",
        f"{task_name}_mRNA_gsea.csv",
        f"{task_name}_survival_analysis.csv",
    ]

    result_files = [
        DownloadableResultFile(
            path=(output_dir / filename).resolve(),
            arcname=f"{archive_root}/{filename}",
            required=False,
        )
        for filename in hybrid_reference_result_filenames
    ]

    deg_filenames = [
        f"{task_name}_{deg_method}_mRNA.csv",
        f"{task_name}_{deg_method}_mRNA_intersect.csv",
        f"{task_name}_{deg_method}_mRNA_venn.csv",
    ]

    for filename in deg_filenames:
        result_files.append(
            DownloadableResultFile(
                path=(output_dir / filename).resolve(),
                arcname=f"{archive_root}/{filename}",
                required=False,
            )
        )

    cm_drug_result_dirname = f"{task_name}_CMdrug_result"

    result_files.append(
        DownloadableResultFile(
            path=(output_dir / cm_drug_result_dirname).resolve(),
            arcname=f"{archive_root}/{cm_drug_result_dirname}/",
            required=False,
            is_directory=True,
        )
    )

    return result_files
