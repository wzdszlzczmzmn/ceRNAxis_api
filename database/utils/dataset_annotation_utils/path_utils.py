from functools import lru_cache
from pathlib import Path
import re
import json

from django.conf import settings

from analysis.utils.workflow_detail_utils.workflow_deg_volcano_utils import WORKFLOW_DEG_SCOPE_ALL, \
    get_workflow_deg_filename

DATASET_NAME_PATTERN = re.compile(r"^[A-Za-z0-9_.-]+$")

RNA_TYPE_SUFFIXES = {
    "mRNA",
    "miRNA",
    "lncRNA",
    "circRNA",
}

TIMEDB_GROUP_BY_SUFFIX_OPTIONS = {
    "grade": {
        "value": "grade",
        "label": "Grade",
    },
    "stage": {
        "value": "stage",
        "label": "Stage",
    },
}

TIMEDB_IGNORED_JSON_GROUP_BY_FIELDS = {
    "c_tumor_grade",
    "c_tumor_stage",
}

TIMEDB_GROUP_TYPES = {
    "other",
    "grade",
    "stage",
}

DATASET_ANNOTATION_AXIS_FINAL_SUFFIX = "_ceRNA_axis_final.csv"

DATASET_ANNOTATION_CERNA_AXIS_SUFFIX = "_ceRNA_axis.csv"

DATASET_ANNOTATION_MAP_IMMUNE_AXIS_SUFFIX = "_map_immune_axis.csv"

DATASET_ANNOTATION_LIMMA_MRNA_SUFFIX = "_limma_mRNA.csv"

DATASET_ANNOTATION_LIMMA_MRNA_INTERSECT_SUFFIX = "_limma_mRNA_intersect.csv"

DATASET_ANNOTATION_CMAP_SUFFIX = "_CMap.csv"

DATASET_ANNOTATION_LOG2FC_BACKGROUND_SUFFIX = "_ceRNA_background.csv"

DATASET_ANNOTATION_EXP_CORRELATION_SUFFIX = "_ceRNA_corr.csv"

DATASET_ANNOTATION_SURVIVAL_SUFFIX = "_survival_analysis.csv"

DATASET_ANNOTATION_MRNA_GSEA_SUFFIX = "_mRNA_gsea.csv"

DATASET_ANNOTATION_SPONGE_RESULT_SUFFIX = "_sponge_result.csv"

DATASET_ANNOTATION_CMDRUG_RESULT_DIR_SUFFIX = "_CMdrug_result"


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


def get_dataset_annotation_sponge_result_file_path(
        *,
        annotation_dir: Path,
        file_prefix: str,
) -> Path:
    return get_dataset_annotation_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
        filename_suffix=DATASET_ANNOTATION_SPONGE_RESULT_SUFFIX,
    )


def get_dataset_annotation_cmdrug_result_dir_path(
        *,
        annotation_dir: Path,
        file_prefix: str,
) -> Path:
    file_prefix = validate_annotation_dataset_name(file_prefix)

    dir_path = (
            annotation_dir / f"{file_prefix}{DATASET_ANNOTATION_CMDRUG_RESULT_DIR_SUFFIX}"
    ).resolve()

    try:
        dir_path.relative_to(annotation_dir)
    except ValueError:
        raise DatasetAnnotationInputError(
            "Invalid dataset annotation CMdrug result directory path."
        )

    return dir_path


def get_dataset_annotation_cerna_axis_file_path(
        *,
        annotation_dir: Path,
        file_prefix: str,
) -> Path:
    return get_dataset_annotation_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
        filename_suffix=DATASET_ANNOTATION_CERNA_AXIS_SUFFIX,
    )


def get_dataset_annotation_map_immune_axis_file_path(
        *,
        annotation_dir: Path,
        file_prefix: str,
) -> Path:
    return get_dataset_annotation_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
        filename_suffix=DATASET_ANNOTATION_MAP_IMMUNE_AXIS_SUFFIX,
    )


def get_dataset_annotation_limma_mrna_file_path(
        *,
        annotation_dir: Path,
        file_prefix: str,
) -> Path:
    return get_dataset_annotation_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
        filename_suffix=DATASET_ANNOTATION_LIMMA_MRNA_SUFFIX,
    )


def get_dataset_annotation_limma_mrna_intersect_file_path(
        *,
        annotation_dir: Path,
        file_prefix: str,
) -> Path:
    return get_dataset_annotation_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
        filename_suffix=DATASET_ANNOTATION_LIMMA_MRNA_INTERSECT_SUFFIX,
    )


def is_existing_file(file_path: Path) -> bool:
    return file_path.exists() and file_path.is_file()


def is_non_empty_file_directory(dir_path: Path) -> bool:
    """
    Return True only when the path is an existing directory containing
    at least one direct child file.

    Nested subdirectories are not searched recursively.
    """
    if not dir_path.exists() or not dir_path.is_dir():
        return False

    try:
        return any(
            child.is_file()
            for child in dir_path.iterdir()
        )
    except OSError:
        return False


def normalize_timedb_json_group_by_label(group_by: str) -> str:
    """
    Convert JSON group_by field to display label.

    Examples:
    - c_tumor_stage -> Tumor Stage
    - c_tumor_grade -> Tumor Grade
    - c_group -> Group
    """
    group_by = str(group_by or "").strip()

    if not group_by:
        return ""

    if group_by.startswith("c_"):
        group_by = group_by[2:]

    return group_by.replace("_", " ").title()


@lru_cache(maxsize=1)
def load_timedb_datasets_info() -> dict:
    json_file = getattr(settings, "TIMEDB_DATASETS_INFO_FILE", None)

    if not json_file:
        raise DatasetAnnotationPathError(
            "TIMEDB_DATASETS_INFO_FILE is not configured."
        )

    json_file = Path(json_file).resolve()

    if not json_file.exists() or not json_file.is_file():
        raise DatasetAnnotationPathError(
            "TIMEDB_datasets_info.json is not available."
        )

    with json_file.open("r", encoding="utf-8") as f:
        raw_data = json.load(f)

    if not isinstance(raw_data, list):
        raise DatasetAnnotationPathError(
            "Invalid TIMEDB_datasets_info.json format."
        )

    result = {}

    for item in raw_data:
        if not isinstance(item, dict):
            continue

        dataset_name = item.get("dataset_name")
        group_by = item.get("group_by", [])

        if not dataset_name:
            continue

        dataset_name = validate_annotation_dataset_name(dataset_name)

        if isinstance(group_by, str):
            group_by = [group_by]

        if not isinstance(group_by, list):
            group_by = []

        result[dataset_name] = [
            str(field).strip()
            for field in group_by
            if str(field or "").strip()
        ]

    return result


def get_timedb_json_group_by_fields(dataset_name: str) -> list[str]:
    dataset_name = validate_annotation_dataset_name(dataset_name)

    datasets_info = load_timedb_datasets_info()

    return datasets_info.get(dataset_name, [])


def get_timedb_other_group_by_option(dataset_name: str) -> dict | None:
    """
    Build the optional TIMEDB ``other`` group from dataset metadata.

    Grade and stage fields are excluded because they are represented by
    dedicated annotation directories. After exclusion, a dataset may contain
    zero or one ``other`` group. More than one field indicates invalid TIMEDB
    metadata and is rejected rather than silently selecting one.
    """
    dataset_name = validate_annotation_dataset_name(dataset_name)

    group_by_fields = []

    for group_by in get_timedb_json_group_by_fields(dataset_name):
        group_by = str(group_by or "").strip()

        if not group_by:
            continue

        if group_by.lower() in TIMEDB_IGNORED_JSON_GROUP_BY_FIELDS:
            continue

        group_by_fields.append(group_by)

    if len(group_by_fields) > 1:
        raise DatasetAnnotationPathError(
            "Invalid TIMEDB dataset metadata: only one other group_by "
            f"field is allowed for dataset {dataset_name}."
        )

    if not group_by_fields:
        return None

    group_by = group_by_fields[0]

    return {
        "value": group_by,
        "label": normalize_timedb_json_group_by_label(group_by),
        "group_type": "other",
        "annotation_dir_name": dataset_name,
        "file_prefix": dataset_name,
    }


def get_timedb_suffix_group_by_option(
        *,
        dataset_name: str,
        suffix: str,
) -> dict | None:
    """
    Suffix annotation directories:
    - module3/{dataset_name}_grade
    - module3/{dataset_name}_stage

    value 直接返回 grade / stage。
    group_type 同样返回 grade / stage，用于后续接口解析目录。
    """
    dataset_name = validate_annotation_dataset_name(dataset_name)
    suffix = str(suffix or "").strip()

    option = TIMEDB_GROUP_BY_SUFFIX_OPTIONS.get(suffix)

    if not option:
        return None

    return {
        "value": option["value"],
        "label": option["label"],
        "group_type": option["value"],
        "source": "folder",
        "annotation_dir_name": f"{dataset_name}_{suffix}",
        "file_prefix": dataset_name,
    }


def build_dataset_annotation_visualization_availability(
        *,
        annotation_dir: Path,
        file_prefix: str,
) -> dict:
    """
    Build visualization availability for one annotation directory.

    ceRNA Annotation Network:
    - {prefix}_ceRNA_axis.csv
    - {prefix}_map_immune_axis.csv

    ceRNA Axis Final Results:
    - {prefix}_ceRNA_axis_final.csv

    CMap Results:
    - {prefix}_CMap.csv

    Expression Volcano Plot:
    - {prefix}_limma_mRNA.csv required
    - {prefix}_limma_mRNA_intersect.csv optional, not used for availability

    Log2FC Correlation Plot:
    - {prefix}_ceRNA_background.csv

    Expression Correlation Plot:
    - {prefix}_ceRNA_corr.csv

    Survival Analysis:
    - {prefix}_survival_analysis.csv

    DEG Pathway Enrichment Plot:
    - {prefix}_mRNA_gsea.csv

    SPONGE:
    - {prefix}_sponge_result.csv

    CMdrug:
    - {prefix}_CMdrug_result/
    - available only when the directory contains at least one file
    """
    cerna_axis_file = get_dataset_annotation_cerna_axis_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
    )

    map_immune_axis_file = get_dataset_annotation_map_immune_axis_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
    )

    axis_final_file = get_dataset_annotation_axis_final_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
    )

    cmap_file = get_dataset_annotation_cmap_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
    )

    limma_mrna_file = get_dataset_annotation_limma_mrna_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
    )

    log2fc_background_file = get_dataset_annotation_log2fc_background_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
    )

    exp_correlation_file = get_dataset_annotation_exp_correlation_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
    )

    survival_file = get_dataset_annotation_survival_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
    )

    mrna_gsea_file = get_dataset_annotation_mrna_gsea_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
    )

    sponge_result_file = get_dataset_annotation_sponge_result_file_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
    )

    cmdrug_result_dir = get_dataset_annotation_cmdrug_result_dir_path(
        annotation_dir=annotation_dir,
        file_prefix=file_prefix,
    )

    visualizations = {
        "annotation_network": (
                is_existing_file(cerna_axis_file)
                and is_existing_file(map_immune_axis_file)
        ),
        "axis_final": is_existing_file(axis_final_file),
        "cmap": is_existing_file(cmap_file),
        "volcano": is_existing_file(limma_mrna_file),
        "log2fc_correlation": is_existing_file(log2fc_background_file),
        "exp_correlation": is_existing_file(exp_correlation_file),
        "survival": is_existing_file(survival_file),
        "deg_pathway": is_existing_file(mrna_gsea_file),
        "sponge": is_existing_file(sponge_result_file),
        "CMdrug": is_non_empty_file_directory(cmdrug_result_dir),
    }

    available_visualization_count = sum(
        1 for available in visualizations.values() if available
    )

    return {
        "visualizations": visualizations,
        "available_visualization_count": available_visualization_count,
    }


def build_timedb_group_by_candidates(dataset_name: str) -> list[dict]:
    """
    Build TIMEDB group candidates in default-selection order:
    other, grade, stage.

    The ``other`` group is optional and unique per dataset.
    """
    dataset_name = validate_annotation_dataset_name(dataset_name)

    candidates = []

    other_option = get_timedb_other_group_by_option(dataset_name)

    if other_option:
        candidates.append(other_option)

    for suffix in ["grade", "stage"]:
        option = get_timedb_suffix_group_by_option(
            dataset_name=dataset_name,
            suffix=suffix,
        )

        if option:
            candidates.append(option)

    return candidates


def validate_timedb_group_selection(
        *,
        dataset_name: str,
        group_by: str,
        group_type: str | None,
) -> dict:
    """
    Validate one TIMEDB Dataset Annotation group selection against
    the same candidate model used by availability.

    Valid combinations are source-dependent:

        other:
            group_by must equal the configured JSON group-by field.

        grade:
            group_by == "grade"

        stage:
            group_by == "stage"

    Returns the matched candidate, including:
        value
        label
        group_type
        annotation_dir_name
        file_prefix
    """
    dataset_name = validate_annotation_dataset_name(
        dataset_name
    )

    group_by = str(
        group_by or ""
    ).strip()

    if not group_by:
        raise DatasetAnnotationInputError(
            "Missing required parameter: group_by."
        )

    group_type = validate_timedb_group_type(
        group_type
    )

    candidates = build_timedb_group_by_candidates(
        dataset_name
    )

    for candidate in candidates:
        if (
                candidate.get("value") == group_by
                and candidate.get("group_type") == group_type
        ):
            return candidate

    allowed = [
        (
            f"{candidate.get('value')}"
            f"/{candidate.get('group_type')}"
        )
        for candidate in candidates
    ]

    allowed_text = (
        ", ".join(allowed)
        if allowed
        else "none"
    )

    raise DatasetAnnotationInputError(
        "Invalid TIMEDB group_by/group_type combination. "
        f"Allowed combinations: {allowed_text}."
    )


def validate_timedb_group_type(group_type: str | None) -> str:
    group_type = str(group_type or "other").strip()

    if not group_type:
        group_type = "other"

    if group_type not in TIMEDB_GROUP_TYPES:
        raise DatasetAnnotationInputError(
            "Invalid group_type parameter. "
            "Allowed values are: other, grade, stage."
        )

    return group_type


def get_timedb_group_type_query(request) -> str:
    return validate_timedb_group_type(
        request.query_params.get("group_type")
    )


def resolve_timedb_group_annotation_dir_name(
        *,
        dataset_name: str,
        group_type: str | None,
) -> str:
    dataset_name = resolve_timedb_annotation_dir_name(dataset_name)
    group_type = validate_timedb_group_type(group_type)

    if group_type == "other":
        return dataset_name

    if group_type == "grade":
        return f"{dataset_name}_grade"

    if group_type == "stage":
        return f"{dataset_name}_stage"

    raise DatasetAnnotationInputError(
        "Invalid group_type parameter."
    )


def resolve_timedb_group_annotation_file_prefix(
        *,
        dataset_name: str,
        group_type: str | None,
) -> str:
    """
    TIMEDB group folders use different directories, but the output file prefix
    remains the original dataset name.

    Examples:
        group_type=other  -> GSE20194/GSE20194_ceRNA_axis.csv
        group_type=grade  -> GSE20194_grade/GSE20194_ceRNA_axis.csv
        group_type=stage  -> GSE20194_stage/GSE20194_ceRNA_axis.csv
    """
    validate_timedb_group_type(group_type)

    return resolve_timedb_annotation_dir_name(dataset_name)
