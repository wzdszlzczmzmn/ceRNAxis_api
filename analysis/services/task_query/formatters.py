from django.utils import timezone

from analysis.utils.paired_cohort_task_utils import (
    get_available_paired_cohort_deg_rna_types,
    get_available_paired_cohort_background_types,
    PairedCohortTaskPathError,
    PairedCohortTaskInputError,
)


def format_datetime(value):
    if value is None:
        return None

    return timezone.localtime(value).strftime("%Y-%m-%d %H:%M:%S")


def normalize_task_rnas_for_response(task) -> dict:
    if not isinstance(task.rnas, dict):
        return {
            "miRNA": [],
            "mRNA": [],
            "lncRNA": [],
            "circRNA": [],
        }

    return {
        "miRNA": task.rnas.get("miRNA", [])
        if isinstance(task.rnas.get("miRNA", []), list)
        else [],
        "mRNA": task.rnas.get("mRNA", [])
        if isinstance(task.rnas.get("mRNA", []), list)
        else [],
        "lncRNA": task.rnas.get("lncRNA", [])
        if isinstance(task.rnas.get("lncRNA", []), list)
        else [],
        "circRNA": task.rnas.get("circRNA", [])
        if isinstance(task.rnas.get("circRNA", []), list)
        else [],
    }


def format_custom_list_query_task(task, position=0) -> dict:
    return {
        "uuid": str(task.uuid),
        "status": task.status,
        "position": position,
        "task_name": task.task_name,
        "map_info": task.map_info,
        "create_time": format_datetime(task.create_time),
        "finish_time": format_datetime(task.finish_time),
        "miRNA_count": task.miRNA_count,
        "mRNA_count": task.mRNA_count,
        "lncRNA_count": task.lncRNA_count,
        "circRNA_count": task.circRNA_count,
        "rnas": normalize_task_rnas_for_response(task),
    }


def get_uploaded_paired_cohort_rna_types(task) -> list[str]:
    """
    Return RNA types whose expression files were actually uploaded/saved.

    mRNA and miRNA are required by the workflow.
    lncRNA and circRNA are optional but at least one should exist.
    """

    uploaded_types = []

    if getattr(task, "mrna_file", ""):
        uploaded_types.append("mRNA")

    if getattr(task, "mirna_file", ""):
        uploaded_types.append("miRNA")

    if getattr(task, "lncrna_file", ""):
        uploaded_types.append("lncRNA")

    if getattr(task, "circrna_file", ""):
        uploaded_types.append("circRNA")

    return uploaded_types


def get_uploaded_paired_cohort_rna_types(task) -> list[str]:
    """
    Return RNA types whose expression files were actually uploaded/saved.

    mRNA and miRNA are required by the workflow.
    lncRNA and circRNA are optional, but at least one of them should exist.
    """

    uploaded_types = []

    if getattr(task, "mrna_file", ""):
        uploaded_types.append("mRNA")

    if getattr(task, "mirna_file", ""):
        uploaded_types.append("miRNA")

    if getattr(task, "lncrna_file", ""):
        uploaded_types.append("lncRNA")

    if getattr(task, "circrna_file", ""):
        uploaded_types.append("circRNA")

    return uploaded_types


def format_paired_cohort_task(task, position=0) -> dict:
    try:
        available_deg_rna_types = get_available_paired_cohort_deg_rna_types(task)
    except (OSError, PairedCohortTaskPathError):
        available_deg_rna_types = []

    try:
        available_background_types = get_available_paired_cohort_background_types(task)
    except (OSError, PairedCohortTaskPathError, PairedCohortTaskInputError):
        available_background_types = []

    uploaded_rna_types = get_uploaded_paired_cohort_rna_types(task)

    return {
        "uuid": str(task.uuid),
        "status": task.status,
        "position": position,
        "task_name": task.task_name,
        "map_info": task.map_info,
        "deg_method": task.deg_method,
        "cancer_type": getattr(task, "cancer_type", ""),
        "use_padj": getattr(task, "use_padj", True),
        "create_time": format_datetime(task.create_time),
        "finish_time": format_datetime(task.finish_time),

        "files": {
            "mrna_file": getattr(task, "mrna_file", ""),
            "mirna_file": getattr(task, "mirna_file", ""),
            "lncrna_file": getattr(task, "lncrna_file", ""),
            "circrna_file": getattr(task, "circrna_file", ""),
            "meta_file": getattr(task, "meta_file", ""),
        },

        "uploaded_rna_types": uploaded_rna_types,
        "has_lncrna_file": bool(getattr(task, "lncrna_file", "")),
        "has_circrna_file": bool(getattr(task, "circrna_file", "")),

        "cutoffs": {
            "mRNA": {
                "logfc_cutoff": task.logfc_cutoff_mrna,
                "pvalue_cutoff": task.padj_cutoff_mrna,
            },
            "miRNA": {
                "logfc_cutoff": task.logfc_cutoff_mirna,
                "pvalue_cutoff": task.padj_cutoff_mirna,
            },
            "lncRNA": {
                "logfc_cutoff": task.logfc_cutoff_lncrna,
                "pvalue_cutoff": task.padj_cutoff_lncrna,
            },
            "circRNA": {
                "logfc_cutoff": task.logfc_cutoff_circrna,
                "pvalue_cutoff": task.padj_cutoff_circrna,
            },
        },

        "available_deg_rna_types": available_deg_rna_types,
        "available_background_types": available_background_types,
    }


def format_hybrid_reference_task(task, position=0) -> dict:
    return {
        "uuid": str(task.uuid),
        "status": task.status,
        "position": position,

        "task_name": task.task_name,
        "map_info": task.map_info,
        "tcga_type": task.tcga_type,
        "lncrna_type": task.lncrna_type,
        "deg_method": task.deg_method,
        "use_padj": getattr(task, "use_padj", True),

        "create_time": format_datetime(task.create_time),
        "finish_time": format_datetime(task.finish_time),

        "files": {
            "mrna_file": getattr(task, "mrna_file", ""),
            "meta_file": getattr(task, "meta_file", ""),
        },

        "uploaded_rna_types": [
            rna_type
            for rna_type, file_value in [
                ("mRNA", getattr(task, "mrna_file", "")),
            ]
            if file_value
        ],

        "cutoffs": {
            "mRNA": {
                "logfc_cutoff": task.logfc_cutoff_mrna,
                "pvalue_cutoff": task.padj_cutoff_mrna,
            },
        },

        "available_deg_rna_types": [
            "mRNA",
        ],
    }
