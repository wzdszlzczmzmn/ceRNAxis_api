from django.utils import timezone

from analysis.utils.paired_cohort_task_utils import (
    get_available_paired_cohort_deg_rna_types,
    get_available_paired_cohort_background_types,
    PairedCohortTaskPathError,
    PairedCohortTaskInputError, PAIRED_COHORT_VALID_BACKGROUND_TYPES,
)
from analysis.utils.workflow_detail_utils.workflow_deg_volcano_utils import get_workflow_available_deg_rna_types, \
    WORKFLOW_DEG_HYBRID_REFERENCE_RNA_TYPES, get_workflow_available_deg_scopes, WORKFLOW_DEG_HYBRID_REFERENCE_SCOPES
from analysis.utils.workflow_detail_utils.workflow_log2fc_background_utils import \
    get_available_workflow_log2fc_background_types, HYBRID_REFERENCE_VALID_BACKGROUND_TYPES


def format_datetime(value):
    if value is None:
        return None

    return timezone.localtime(value).strftime("%Y-%m-%d %H:%M:%S")


CUSTOM_LIST_QUERY_RNA_TYPES = (
    "miRNA",
    "mRNA",
    "mRNA_up",
    "mRNA_down",
    "lncRNA",
    "circRNA",
)


def normalize_task_rnas_for_response(task) -> dict:
    """
    Normalize the CustomListQueryTask.rnas JSON value.

    New tasks support:
    - Non-directional mRNA input: mRNA
    - Directional mRNA input: mRNA_up and mRNA_down

    Missing keys in legacy tasks are returned as empty lists.
    Invalid non-list values are also normalized to empty lists.
    """
    raw_rnas = task.rnas

    if not isinstance(raw_rnas, dict):
        raw_rnas = {}

    return {
        rna_type: (
            raw_rnas.get(rna_type, [])
            if isinstance(raw_rnas.get(rna_type, []), list)
            else []
        )
        for rna_type in CUSTOM_LIST_QUERY_RNA_TYPES
    }


def format_custom_list_query_task(task, position=0) -> dict:
    rnas = normalize_task_rnas_for_response(task)

    miRNA_count = len(rnas["miRNA"])
    mRNA_count = len(rnas["mRNA"])
    mRNA_up_count = len(rnas["mRNA_up"])
    mRNA_down_count = len(rnas["mRNA_down"])
    lncRNA_count = len(rnas["lncRNA"])
    circRNA_count = len(rnas["circRNA"])

    total_rna_count = (
        miRNA_count
        + mRNA_count
        + mRNA_up_count
        + mRNA_down_count
        + lncRNA_count
        + circRNA_count
    )

    has_mrna_direction = bool(
        getattr(task, "has_mrna_direction", False)
    )

    return {
        "uuid": str(task.uuid),
        "status": task.status,
        "status_label": task.get_status_display(),
        "position": position,

        "task_name": task.task_name,
        "user": getattr(task, "user", "") or "",

        # Deprecated for CustomListQueryTask, retained for compatibility.
        "map_info": getattr(task, "map_info", "") or "",

        # New Module 1 fields.
        "cancer_type": getattr(task, "cancer_type", "") or "",
        "has_mrna_direction": has_mrna_direction,

        # Explicitly describes which mRNA fields the task uses.
        "mrna_input_mode": (
            "directional"
            if has_mrna_direction
            else "standard"
        ),

        "create_time": format_datetime(task.create_time),
        "finish_time": format_datetime(task.finish_time),

        # Retained top-level count fields for backward compatibility.
        "miRNA_count": miRNA_count,
        "mRNA_count": mRNA_count,
        "mRNA_up_count": mRNA_up_count,
        "mRNA_down_count": mRNA_down_count,
        "lncRNA_count": lncRNA_count,
        "circRNA_count": circRNA_count,
        "total_rna_count": total_rna_count,

        "rna_counts": {
            "miRNA": miRNA_count,
            "mRNA": mRNA_count,
            "mRNA_up": mRNA_up_count,
            "mRNA_down": mRNA_down_count,
            "lncRNA": lncRNA_count,
            "circRNA": circRNA_count,
            "total": total_rna_count,
        },

        "rnas": rnas,
    }


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

    available_background_types = get_available_workflow_log2fc_background_types(
        task=task,
        valid_types=PAIRED_COHORT_VALID_BACKGROUND_TYPES,
    )

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
    available_background_types = get_available_workflow_log2fc_background_types(
        task=task,
        valid_types=HYBRID_REFERENCE_VALID_BACKGROUND_TYPES,
    )

    available_deg_rna_types = get_workflow_available_deg_rna_types(
        task=task,
        valid_rna_types=WORKFLOW_DEG_HYBRID_REFERENCE_RNA_TYPES,
    )

    available_deg_scopes = get_workflow_available_deg_scopes(
        task=task,
        rna_type="mRNA",
        valid_scopes=WORKFLOW_DEG_HYBRID_REFERENCE_SCOPES,
    )

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

        "available_deg_rna_types": available_deg_rna_types,
        "available_deg_scopes": available_deg_scopes,
        "available_background_types": available_background_types,
    }
