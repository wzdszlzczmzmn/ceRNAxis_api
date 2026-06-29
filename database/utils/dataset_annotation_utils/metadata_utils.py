from pathlib import Path

from analysis.utils.workflow_detail_utils.workflow_deg_volcano_utils import (
    WORKFLOW_DEG_SCOPE_ALL,
    WORKFLOW_DEG_HYBRID_REFERENCE_RNA_TYPES,
    WORKFLOW_DEG_HYBRID_REFERENCE_SCOPES,
)

from analysis.utils.workflow_detail_utils.workflow_log2fc_background_utils import (
    HYBRID_REFERENCE_VALID_BACKGROUND_TYPES,
    PAIRED_COHORT_VALID_BACKGROUND_TYPES,
    get_workflow_available_background_types,
    read_log2fc_background_file_by_path,
    WorkflowLog2FCBackgroundInputError,
)

from analysis.utils.workflow_detail_utils.workflow_deg_volcano_utils import (
    WORKFLOW_DEG_PAIRED_COHORT_RNA_TYPES,
    WORKFLOW_DEG_PAIRED_COHORT_SCOPES,
)

from database.utils.dataset_annotation_utils.path_utils import (
    get_dataset_annotation_deg_file_path,
    get_dataset_annotation_log2fc_background_file_path,
)


DEFAULT_DEG_METHOD = "limma"
DEFAULT_USE_PADJ = False


DEFAULT_TCGA_CUTOFFS = {
    "mRNA": {
        "logfc_cutoff": 0.000001,
        "pvalue_cutoff": 0.05,
    },
    "miRNA": {
        "logfc_cutoff": 0.000001,
        "pvalue_cutoff": 0.05,
    },
    "lncRNA": {
        "logfc_cutoff": 0.000001,
        "pvalue_cutoff": 0.05,
    },
    "circRNA": {
        "logfc_cutoff": 0.000001,
        "pvalue_cutoff": 0.05,
    },
}


DEFAULT_TIMEDB_CUTOFFS = {
    "mRNA": {
        "logfc_cutoff": 0.000001,
        "pvalue_cutoff": 0.05,
    },
}


def get_dataset_annotation_available_deg_rna_types(
    *,
    annotation_dir: Path,
    file_prefix: str,
    deg_method: str,
    valid_rna_types: list[str],
    deg_scope: str = WORKFLOW_DEG_SCOPE_ALL,
) -> list[str]:
    available_rna_types = []

    for rna_type in valid_rna_types:
        file_path = get_dataset_annotation_deg_file_path(
            annotation_dir=annotation_dir,
            file_prefix=file_prefix,
            deg_method=deg_method,
            rna_type=rna_type,
            deg_scope=deg_scope,
        )

        if file_path.exists() and file_path.is_file():
            available_rna_types.append(rna_type)

    return available_rna_types


def get_dataset_annotation_available_deg_scopes(
    *,
    annotation_dir: Path,
    file_prefix: str,
    deg_method: str,
    rna_type: str,
    valid_scopes: list[str],
) -> list[str]:
    available_scopes = []

    for deg_scope in valid_scopes:
        file_path = get_dataset_annotation_deg_file_path(
            annotation_dir=annotation_dir,
            file_prefix=file_prefix,
            deg_method=deg_method,
            rna_type=rna_type,
            deg_scope=deg_scope,
        )

        if file_path.exists() and file_path.is_file():
            available_scopes.append(deg_scope)

    return available_scopes


def get_dataset_annotation_available_background_types(
    *,
    annotation_dir: Path,
    file_prefix: str,
    valid_types: list[str],
) -> list[str]:
    background_file = get_dataset_annotation_log2fc_background_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
    )

    try:
        _, df = read_log2fc_background_file_by_path(background_file)
    except (
        FileNotFoundError,
        WorkflowLog2FCBackgroundInputError,
    ):
        return []

    return get_workflow_available_background_types(
        df=df,
        valid_types=valid_types,
    )


def build_tcga_dataset_annotation_metadata(
    *,
    annotation_dir: Path,
    file_prefix: str,
    deg_method: str = DEFAULT_DEG_METHOD,
) -> dict:
    """
    TCGA Annotation comes from Paired Cohort.
    """

    available_deg_rna_types = get_dataset_annotation_available_deg_rna_types(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
        deg_method=deg_method,
        valid_rna_types=WORKFLOW_DEG_PAIRED_COHORT_RNA_TYPES,
        deg_scope=WORKFLOW_DEG_SCOPE_ALL,
    )

    available_background_types = get_dataset_annotation_available_background_types(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
        valid_types=PAIRED_COHORT_VALID_BACKGROUND_TYPES,
    )

    return {
        "network_source_task_type": "PairedCohortTask",
        "deg_method": deg_method,
        "use_padj": DEFAULT_USE_PADJ,
        "cutoffs": DEFAULT_TCGA_CUTOFFS,
        "available_deg_rna_types": available_deg_rna_types,
        "available_deg_scopes": WORKFLOW_DEG_PAIRED_COHORT_SCOPES,
        "available_background_types": available_background_types,
    }


def build_timedb_dataset_annotation_metadata(
    *,
    annotation_dir: Path,
    file_prefix: str,
    deg_method: str = DEFAULT_DEG_METHOD,
) -> dict:
    """
    TIMEDB Annotation comes from Module3 / Hybrid Reference.

    Module3 DEG semantics:
        valid RNA types: mRNA
        valid scopes: all, intersect
    """

    available_deg_rna_types = get_dataset_annotation_available_deg_rna_types(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
        deg_method=deg_method,
        valid_rna_types=WORKFLOW_DEG_HYBRID_REFERENCE_RNA_TYPES,
        deg_scope=WORKFLOW_DEG_SCOPE_ALL,
    )

    default_rna_type = (
        "mRNA"
        if "mRNA" in available_deg_rna_types
        else available_deg_rna_types[0]
        if available_deg_rna_types
        else None
    )

    available_deg_scopes = (
        get_dataset_annotation_available_deg_scopes(
            annotation_dir=annotation_dir,
            file_prefix=file_prefix,
            deg_method=deg_method,
            rna_type=default_rna_type,
            valid_scopes=WORKFLOW_DEG_HYBRID_REFERENCE_SCOPES,
        )
        if default_rna_type
        else []
    )

    available_background_types = get_dataset_annotation_available_background_types(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
        valid_types=HYBRID_REFERENCE_VALID_BACKGROUND_TYPES,
    )

    return {
        "network_source_task_type": "HybridReferenceTask",
        "deg_method": deg_method,
        "use_padj": DEFAULT_USE_PADJ,
        "cutoffs": DEFAULT_TIMEDB_CUTOFFS,
        "available_deg_rna_types": available_deg_rna_types,
        "available_deg_scopes": available_deg_scopes,
        "available_background_types": available_background_types,
    }


def build_dataset_annotation_metadata(
    *,
    annotation_dir: Path,
    file_prefix: str,
    source: str,
    deg_method: str = DEFAULT_DEG_METHOD,
) -> dict:
    if source == "TCGA":
        return build_tcga_dataset_annotation_metadata(
            annotation_dir=annotation_dir,
            file_prefix=file_prefix,
            deg_method=deg_method,
        )

    if source == "TIMEDB":
        return build_timedb_dataset_annotation_metadata(
            annotation_dir=annotation_dir,
            file_prefix=file_prefix,
            deg_method=deg_method,
        )

    raise ValueError(f"Unsupported annotation source: {source}")
