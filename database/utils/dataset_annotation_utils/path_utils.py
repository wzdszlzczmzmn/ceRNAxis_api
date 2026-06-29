from pathlib import Path
import re

from analysis.utils.workflow_detail_utils.workflow_deg_volcano_utils import WORKFLOW_DEG_SCOPE_ALL, \
    get_workflow_deg_filename

DATASET_NAME_PATTERN = re.compile(r"^[A-Za-z0-9_.-]+$")

RNA_TYPE_SUFFIXES = {
    "mRNA",
    "miRNA",
    "lncRNA",
    "circRNA",
}


DATASET_ANNOTATION_AXIS_FINAL_SUFFIX = "_ceRNA_axis_final.csv"


DATASET_ANNOTATION_CMAP_SUFFIX = "_CMap.csv"


DATASET_ANNOTATION_LOG2FC_BACKGROUND_SUFFIX = "_ceRNA_background.csv"


DATASET_ANNOTATION_EXP_CORRELATION_SUFFIX = "_ceRNA_corr.csv"


DATASET_ANNOTATION_SURVIVAL_SUFFIX = "_survival_analysis.csv"


DATASET_ANNOTATION_MRNA_GSEA_SUFFIX = "_mRNA_gsea.csv"


class DatasetAnnotationInputError(ValueError):
    pass


class DatasetAnnotationPathError(ValueError):
    pass


def validate_annotation_dataset_name(dataset_name: str) -> str:
    dataset_name = str(dataset_name or "").strip()

    if not dataset_name:
        raise DatasetAnnotationInputError(
            "Missing required parameter: dataset."
        )

    if "/" in dataset_name or "\\" in dataset_name or ".." in dataset_name:
        raise DatasetAnnotationInputError(
            "Invalid dataset parameter."
        )

    if dataset_name.startswith("."):
        raise DatasetAnnotationInputError(
            "Invalid dataset parameter."
        )

    if not DATASET_NAME_PATTERN.match(dataset_name):
        raise DatasetAnnotationInputError(
            "Invalid dataset parameter. "
            "Only letters, numbers, underscore, hyphen and dot are allowed."
        )

    return dataset_name


def strip_rna_type_suffix(dataset_name: str) -> str:
    dataset_name = validate_annotation_dataset_name(dataset_name)

    parts = dataset_name.rsplit("_", 1)

    if len(parts) != 2:
        return dataset_name

    prefix, suffix = parts

    if suffix in RNA_TYPE_SUFFIXES:
        return prefix

    return dataset_name


def resolve_tcga_annotation_dir_name(dataset_name: str) -> str:
    """
    TCGA_ACC_mRNA -> TCGA_ACC
    """
    return strip_rna_type_suffix(dataset_name)


def resolve_timedb_annotation_dir_name(dataset_name: str) -> str:
    """
    目前保持原名。

    如果你的 TIMEDB 目录也是 GSE19750 而不是 GSE19750_mRNA，
    就改成 return strip_rna_type_suffix(dataset_name)
    """
    return validate_annotation_dataset_name(dataset_name)


def get_dataset_query_name(request) -> str:
    dataset_name = request.query_params.get("dataset")

    if dataset_name is None:
        dataset_name = request.query_params.get("dataset_name")

    return validate_annotation_dataset_name(dataset_name)


def resolve_dataset_annotation_dir(
    *,
    annotation_root_dir,
    annotation_dir_name: str,
) -> Path:
    annotation_dir_name = validate_annotation_dataset_name(annotation_dir_name)

    root_dir = Path(annotation_root_dir).resolve()

    if not root_dir.exists() or not root_dir.is_dir():
        raise DatasetAnnotationPathError(
            "Dataset annotation root directory is not available."
        )

    annotation_dir = (root_dir / annotation_dir_name).resolve()

    try:
        annotation_dir.relative_to(root_dir)
    except ValueError:
        raise DatasetAnnotationInputError(
            "Invalid dataset annotation directory."
        )

    return annotation_dir


def get_dataset_annotation_file_path(
    *,
    annotation_dir: Path,
    file_prefix: str,
    filename_suffix: str,
) -> Path:
    file_prefix = validate_annotation_dataset_name(file_prefix)

    file_path = (
        annotation_dir / f"{file_prefix}{filename_suffix}"
    ).resolve()

    try:
        file_path.relative_to(annotation_dir)
    except ValueError:
        raise DatasetAnnotationInputError(
            "Invalid dataset annotation file path."
        )

    return file_path


def get_dataset_annotation_axis_final_file_path(
    *,
    annotation_dir: Path,
    file_prefix: str,
) -> Path:
    return get_dataset_annotation_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
        filename_suffix=DATASET_ANNOTATION_AXIS_FINAL_SUFFIX,
    )


def get_dataset_annotation_cmap_file_path(
    *,
    annotation_dir: Path,
    file_prefix: str,
) -> Path:
    return get_dataset_annotation_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
        filename_suffix=DATASET_ANNOTATION_CMAP_SUFFIX,
    )


def get_dataset_annotation_deg_file_path(
    *,
    annotation_dir: Path,
    file_prefix: str,
    deg_method: str,
    rna_type: str,
    deg_scope: str = WORKFLOW_DEG_SCOPE_ALL,
) -> Path:
    filename = get_workflow_deg_filename(
        task_name=file_prefix,
        deg_method=deg_method,
        rna_type=rna_type,
        deg_scope=deg_scope,
    )

    file_path = (annotation_dir / filename).resolve()

    try:
        file_path.relative_to(annotation_dir)
    except ValueError:
        raise DatasetAnnotationInputError(
            "Invalid dataset annotation DEG file path."
        )

    return file_path


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
        try:
            file_path = get_dataset_annotation_deg_file_path(
                annotation_dir=annotation_dir,
                file_prefix=file_prefix,
                deg_method=deg_method,
                rna_type=rna_type,
                deg_scope=deg_scope,
            )
        except Exception:
            continue

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
        try:
            file_path = get_dataset_annotation_deg_file_path(
                annotation_dir=annotation_dir,
                file_prefix=file_prefix,
                deg_method=deg_method,
                rna_type=rna_type,
                deg_scope=deg_scope,
            )
        except Exception:
            continue

        if file_path.exists() and file_path.is_file():
            available_scopes.append(deg_scope)

    return available_scopes


def get_dataset_annotation_log2fc_background_file_path(
    *,
    annotation_dir: Path,
    file_prefix: str,
) -> Path:
    return get_dataset_annotation_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
        filename_suffix=DATASET_ANNOTATION_LOG2FC_BACKGROUND_SUFFIX,
    )


def get_dataset_annotation_exp_correlation_file_path(
    *,
    annotation_dir: Path,
    file_prefix: str,
) -> Path:
    return get_dataset_annotation_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
        filename_suffix=DATASET_ANNOTATION_EXP_CORRELATION_SUFFIX,
    )


def get_dataset_annotation_survival_file_path(
    *,
    annotation_dir: Path,
    file_prefix: str,
) -> Path:
    return get_dataset_annotation_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
        filename_suffix=DATASET_ANNOTATION_SURVIVAL_SUFFIX,
    )


def get_dataset_annotation_mrna_gsea_file_path(
    *,
    annotation_dir: Path,
    file_prefix: str,
) -> Path:
    return get_dataset_annotation_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
        filename_suffix=DATASET_ANNOTATION_MRNA_GSEA_SUFFIX,
    )
