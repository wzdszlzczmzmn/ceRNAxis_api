from functools import lru_cache
from pathlib import Path
import json
import re

from django.conf import settings

from database.utils.dataset_annotation_utils.path_utils import (
    DatasetAnnotationInputError,
    DatasetAnnotationPathError,
    resolve_dataset_annotation_dir,
    validate_annotation_dataset_name,
)


SCST_DATA_TYPES = {
    "sc",
    "st",
}


SCST_DATASET_ID_COLUMN_MAP = {
    "sc": "cell_id",
    "st": "spot_id",
}


SCST_DATASET_ANNOTATION_ROOT_SETTING_MAP = {
    "sc": "TISCH2_DATASET_ANNOTATIONS_DIR",
    "st": "SCTML_DATASET_ANNOTATIONS_DIR",
}


SCST_DATASET_META_ROOT_SETTING_MAP = {
    "sc": "TISCH2_DATASET_BASE_DIR",
    "st": "SCTML_DATASET_BASE_DIR",
}


SCST_DATASET_GROUP_COL_SETTING_NAME = (
    "SCST_DATASET_GROUP_COL_FILE"
)

SCST_DATASET_META_FILENAME_SUFFIX = "_meta.csv"


SCST_GROUP_BY_META_COLUMN_MAP = {
    "celltype malignancy": "Celltype (malignancy)",
    "Celltype major lineage": "Celltype (major-lineage)",
    "Celltype minor lineage": "Celltype (minor-lineage)",
    "celllabel": "c_cell_label",
    "celltype": "c_cell_type",
}


def validate_scst_data_type(
    data_type: str | None,
) -> str:
    data_type = str(data_type or "").strip().lower()

    if not data_type:
        raise DatasetAnnotationInputError(
            "Missing required parameter: data_type."
        )

    if data_type not in SCST_DATA_TYPES:
        raise DatasetAnnotationInputError(
            "Invalid data_type parameter. "
            "Allowed values are: sc, st."
        )

    return data_type


def get_scst_data_type_query(request) -> str:
    return validate_scst_data_type(
        request.query_params.get("data_type")
    )


def get_scst_dataset_id_column(
    data_type: str,
) -> str:
    data_type = validate_scst_data_type(
        data_type
    )

    return SCST_DATASET_ID_COLUMN_MAP[
        data_type
    ]


def get_required_scst_setting_path(
    setting_name: str,
) -> Path:
    value = getattr(
        settings,
        setting_name,
        None,
    )

    if not value:
        raise DatasetAnnotationPathError(
            f"{setting_name} is not configured."
        )

    return Path(value).resolve()


def get_scst_dataset_annotation_root_dir(
    data_type: str,
) -> Path:
    """
    Return the annotation-result root for one SC/ST source.

    SC:
        settings.TISCH2_DATASET_ANNOTATIONS_DIR

    ST:
        settings.SCTML_DATASET_ANNOTATIONS_DIR
    """
    data_type = validate_scst_data_type(
        data_type
    )

    setting_name = (
        SCST_DATASET_ANNOTATION_ROOT_SETTING_MAP[
            data_type
        ]
    )

    root_dir = get_required_scst_setting_path(
        setting_name
    )

    if not root_dir.exists() or not root_dir.is_dir():
        raise DatasetAnnotationPathError(
            f"{setting_name} is not available."
        )

    return root_dir


def get_scst_dataset_meta_root_dir(
    data_type: str,
) -> Path:
    """
    Return the dataset metadata root for one SC/ST source.

    SC:
        settings.TISCH2_DATASET_BASE_DIR

    ST:
        settings.SCTML_DATASET_BASE_DIR
    """
    data_type = validate_scst_data_type(
        data_type
    )

    setting_name = (
        SCST_DATASET_META_ROOT_SETTING_MAP[
            data_type
        ]
    )

    root_dir = get_required_scst_setting_path(
        setting_name
    )

    if not root_dir.exists() or not root_dir.is_dir():
        raise DatasetAnnotationPathError(
            f"{setting_name} is not available."
        )

    return root_dir


def get_scst_dataset_group_col_file_path() -> Path:
    json_file = get_required_scst_setting_path(
        SCST_DATASET_GROUP_COL_SETTING_NAME
    )

    if not json_file.exists() or not json_file.is_file():
        raise DatasetAnnotationPathError(
            "SCST_DATASET_GROUP_COL_FILE is not available."
        )

    return json_file


@lru_cache(maxsize=1)
def load_scst_dataset_group_cols() -> dict[str, list[str]]:
    """
    Load:
        dataset_name -> list[group_by column]

    from:
        settings.SCST_DATASET_GROUP_COL_FILE
    """
    json_file = (
        get_scst_dataset_group_col_file_path()
    )

    try:
        with json_file.open(
            "r",
            encoding="utf-8",
        ) as file_obj:
            raw_data = json.load(file_obj)

    except (
        UnicodeDecodeError,
        json.JSONDecodeError,
        OSError,
    ) as exc:
        raise DatasetAnnotationPathError(
            "Failed to read SC/ST dataset group-by "
            "configuration."
        ) from exc

    if not isinstance(raw_data, dict):
        raise DatasetAnnotationPathError(
            "Invalid SC/ST dataset group-by "
            "configuration format."
        )

    result: dict[str, list[str]] = {}

    for (
        raw_dataset_name,
        raw_group_by_fields,
    ) in raw_data.items():
        dataset_name = (
            validate_annotation_dataset_name(
                raw_dataset_name
            )
        )

        if isinstance(
            raw_group_by_fields,
            str,
        ):
            raw_group_by_fields = [
                raw_group_by_fields
            ]

        if not isinstance(
            raw_group_by_fields,
            list,
        ):
            raise DatasetAnnotationPathError(
                "Invalid SC/ST group-by configuration "
                f"for dataset {dataset_name}."
            )

        group_by_fields = []
        seen = set()

        for raw_group_by in raw_group_by_fields:
            group_by = str(
                raw_group_by or ""
            ).strip()

            if not group_by:
                continue

            if group_by in seen:
                continue

            group_by_fields.append(
                group_by
            )
            seen.add(group_by)

        result[dataset_name] = (
            group_by_fields
        )

    return result


def get_scst_dataset_group_by_fields(
    dataset_name: str,
) -> list[str]:
    dataset_name = (
        validate_annotation_dataset_name(
            dataset_name
        )
    )

    return list(
        load_scst_dataset_group_cols().get(
            dataset_name,
            [],
        )
    )


def get_scst_group_by_query(
    request,
) -> str:
    """
    Return the required public/configured SC/ST group_by value.

    Dataset-specific membership validation is intentionally kept
    separate because it requires the resolved dataset name.
    """
    group_by = str(
        request.query_params.get(
            "group_by",
            "",
        )
        or ""
    ).strip()

    if not group_by:
        raise DatasetAnnotationInputError(
            "Missing required parameter: group_by."
        )

    return group_by


def validate_scst_dataset_group_by(
    *,
    dataset_name: str,
    group_by: str,
) -> str:
    """
    Validate that group_by is configured for the requested dataset.

    Matching is exact because configured group_by values are public
    API/business values and are also used to resolve result folders.
    """
    dataset_name = (
        validate_annotation_dataset_name(
            dataset_name
        )
    )

    group_by = str(
        group_by or ""
    ).strip()

    if not group_by:
        raise DatasetAnnotationInputError(
            "Missing required parameter: group_by."
        )

    configured_group_by_fields = (
        get_scst_dataset_group_by_fields(
            dataset_name
        )
    )

    if group_by not in configured_group_by_fields:
        raise DatasetAnnotationInputError(
            "Invalid SC/ST group_by for this dataset."
        )

    return group_by


def get_scst_group_by_meta_column(
    group_by: str,
) -> str:
    """
    Resolve a configured SC/ST group-by name to the actual
    column name used in the dataset metadata CSV.

    The configured group-by name remains the public/API value
    and is also used to build the annotation result directory.

    Examples:
        "celltype malignancy"
            -> "Celltype (malignancy)"

        "Celltype major lineage"
            -> "Celltype (major-lineage)"
    """
    group_by = str(
        group_by or ""
    ).strip()

    if not group_by:
        raise DatasetAnnotationInputError(
            "Missing SC/ST group_by value."
        )

    try:
        return SCST_GROUP_BY_META_COLUMN_MAP[
            group_by
        ]
    except KeyError as exc:
        raise DatasetAnnotationPathError(
            "Unsupported SC/ST configured group_by: "
            f"{group_by}."
        ) from exc


def normalize_scst_group_by_dir_suffix(
    group_by: str,
) -> str:
    """
    Convert a group-by column into the result-directory
    suffix.

    Example:
        Celltype major lineage
            -> Celltype_major_lineage

    Only whitespace normalization is performed.
    """
    group_by = str(
        group_by or ""
    ).strip()

    if not group_by:
        raise DatasetAnnotationInputError(
            "Missing SC/ST group_by value."
        )

    return re.sub(
        r"\s+",
        "_",
        group_by,
    )


def get_scst_dataset_group_annotation_dir_name(
    *,
    dataset_name: str,
    group_by: str,
) -> str:
    dataset_name = (
        validate_annotation_dataset_name(
            dataset_name
        )
    )

    group_by_suffix = (
        normalize_scst_group_by_dir_suffix(
            group_by
        )
    )

    annotation_dir_name = (
        f"{dataset_name}_{group_by_suffix}"
    )

    # Reuse the generic dataset/path-safe name validator.
    return validate_annotation_dataset_name(
        annotation_dir_name
    )


def resolve_scst_dataset_group_annotation_dir(
    *,
    dataset_name: str,
    group_by: str,
    data_type: str,
) -> Path:
    """
    Resolve one group-by annotation-result directory.

    SC currently resolves under:
        TISCH2_DATASET_ANNOTATIONS_DIR/
            {dataset_name}_{normalized_group_by}/

    ST is intentionally kept as its own branch so its
    storage rule can diverge later without affecting SC.
    """
    data_type = validate_scst_data_type(
        data_type
    )

    annotation_dir_name = (
        get_scst_dataset_group_annotation_dir_name(
            dataset_name=dataset_name,
            group_by=group_by,
        )
    )

    if data_type == "sc":
        annotation_root_dir = (
            get_scst_dataset_annotation_root_dir(
                "sc"
            )
        )

        return resolve_dataset_annotation_dir(
            annotation_root_dir=annotation_root_dir,
            annotation_dir_name=annotation_dir_name,
        )

    if data_type == "st":
        annotation_root_dir = (
            get_scst_dataset_annotation_root_dir(
                "st"
            )
        )

        # Reserved ST branch:
        # change only this block if SCTML annotation
        # directories later use a different layout.
        return resolve_dataset_annotation_dir(
            annotation_root_dir=annotation_root_dir,
            annotation_dir_name=annotation_dir_name,
        )

    raise DatasetAnnotationInputError(
        "Invalid SC/ST data_type."
    )


def get_scst_dataset_meta_file_path(
    *,
    dataset_name: str,
    data_type: str,
) -> Path:
    """
    Resolve the dataset-level metadata CSV.

    SC:
        TISCH2_DATASET_BASE_DIR/
            {dataset_name}_meta.csv

    ST:
        SCTML_DATASET_BASE_DIR/
            {dataset_name}_meta.csv
    """
    data_type = validate_scst_data_type(
        data_type
    )

    dataset_name = (
        validate_annotation_dataset_name(
            dataset_name
        )
    )

    if data_type == "sc":
        meta_root_dir = (
            get_scst_dataset_meta_root_dir(
                "sc"
            )
        )

    elif data_type == "st":
        meta_root_dir = (
            get_scst_dataset_meta_root_dir(
                "st"
            )
        )

    else:
        raise DatasetAnnotationInputError(
            "Invalid SC/ST data_type."
        )

    meta_file = (
        meta_root_dir
        / (
            f"{dataset_name}"
            f"{SCST_DATASET_META_FILENAME_SUFFIX}"
        )
    ).resolve()

    try:
        meta_file.relative_to(
            meta_root_dir
        )
    except ValueError as exc:
        raise DatasetAnnotationInputError(
            "Invalid SC/ST dataset metadata file path."
        ) from exc

    return meta_file

def should_skip_scst_group_value_for_results(
    group_value: str,
) -> bool:
    """
    Return True when a metadata group value should be ignored
    for Dataset Annotation result availability.

    Current rule:
        values containing "/" are skipped because the upstream
        SC/ST annotation workflow does not generate result files
        for those values.

    The original metadata row still contributes to sample_count.
    """
    group_value = str(
        group_value or ""
    ).strip()

    return "/" in group_value


def get_scst_group_value_query(
    request,
) -> str:
    """
    Return the required SC/ST group_value used to resolve one
    visualization result.

    Values containing "/" are intentionally unavailable because
    the upstream annotation workflow does not generate result
    files for those metadata values.
    """
    group_value = str(
        request.query_params.get(
            "group_value",
            "",
        )
        or ""
    ).strip()

    if not group_value:
        raise DatasetAnnotationInputError(
            "Missing required parameter: group_value."
        )

    if should_skip_scst_group_value_for_results(
        group_value
    ):
        raise DatasetAnnotationInputError(
            "SC/ST group_value is unavailable for result files."
        )

    return validate_scst_group_value_path_component(
        group_value
    )


def validate_scst_group_value_path_component(
    group_value: str,
) -> str:
    """
    Validate a group value before embedding it in an SC/ST
    annotation result filename.

    Group values are intentionally not normalized because the
    result files use the original metadata value verbatim.

    Slash-containing values are filtered before this function.
    Backslashes and traversal-like sequences are rejected here.
    """
    group_value = str(
        group_value or ""
    ).strip()

    if not group_value:
        raise DatasetAnnotationInputError(
            "Missing SC/ST group_value."
        )

    if (
        "\\" in group_value
        or ".." in group_value
    ):
        raise DatasetAnnotationInputError(
            "Invalid SC/ST group_value for result path."
        )

    return group_value


def get_scst_group_value_result_path(
    *,
    annotation_dir: Path,
    filename: str,
) -> Path:
    """
    Resolve a direct child result path and enforce that it
    remains inside the group-by annotation directory.
    """
    annotation_dir = Path(
        annotation_dir
    ).resolve()

    result_path = (
        annotation_dir / filename
    ).resolve()

    try:
        result_path.relative_to(
            annotation_dir
        )
    except ValueError as exc:
        raise DatasetAnnotationInputError(
            "Invalid SC/ST annotation result path."
        ) from exc

    return result_path


def get_scst_dataset_cerna_axis_file_path(
    *,
    annotation_dir: Path,
    dataset_name: str,
    group_value: str,
) -> Path:
    dataset_name = (
        validate_annotation_dataset_name(
            dataset_name
        )
    )
    group_value = (
        validate_scst_group_value_path_component(
            group_value
        )
    )

    return get_scst_group_value_result_path(
        annotation_dir=annotation_dir,
        filename=(
            f"{dataset_name}_ceRNA_"
            f"{group_value}_axis.csv"
        ),
    )


def get_scst_dataset_map_immune_axis_file_path(
    *,
    annotation_dir: Path,
    dataset_name: str,
    group_value: str,
) -> Path:
    dataset_name = (
        validate_annotation_dataset_name(
            dataset_name
        )
    )
    group_value = (
        validate_scst_group_value_path_component(
            group_value
        )
    )

    return get_scst_group_value_result_path(
        annotation_dir=annotation_dir,
        filename=(
            f"{dataset_name}_map_immune_axis_"
            f"{group_value}.csv"
        ),
    )


def get_scst_dataset_axis_final_file_path(
    *,
    annotation_dir: Path,
    dataset_name: str,
    group_value: str,
) -> Path:
    dataset_name = (
        validate_annotation_dataset_name(
            dataset_name
        )
    )
    group_value = (
        validate_scst_group_value_path_component(
            group_value
        )
    )

    return get_scst_group_value_result_path(
        annotation_dir=annotation_dir,
        filename=(
            f"{dataset_name}_ceRNA_"
            f"{group_value}_axis_final.csv"
        ),
    )


def get_scst_dataset_cmap_file_path(
    *,
    annotation_dir: Path,
    dataset_name: str,
    group_value: str,
) -> Path:
    dataset_name = (
        validate_annotation_dataset_name(
            dataset_name
        )
    )
    group_value = (
        validate_scst_group_value_path_component(
            group_value
        )
    )

    return get_scst_group_value_result_path(
        annotation_dir=annotation_dir,
        filename=(
            f"{dataset_name}_CMap_"
            f"{group_value}.csv"
        ),
    )


def get_scst_dataset_deg_file_path(
    *,
    annotation_dir: Path,
    dataset_name: str,
    group_value: str,
) -> Path:
    """
    SC/ST Dataset Annotation currently has one DEG file per
    group value and supports mRNA only.

    Expected filename:
        {dataset_name}_deg_{group_value}.csv
    """
    dataset_name = (
        validate_annotation_dataset_name(
            dataset_name
        )
    )
    group_value = (
        validate_scst_group_value_path_component(
            group_value
        )
    )

    return get_scst_group_value_result_path(
        annotation_dir=annotation_dir,
        filename=(
            f"{dataset_name}_deg_"
            f"{group_value}.csv"
        ),
    )


def get_scst_dataset_deg_intersect_file_path(
    *,
    annotation_dir: Path,
    dataset_name: str,
    group_value: str,
) -> Path:
    """
    Resolve the SC/ST Dataset Annotation intersect DEG file.

    Expected filename:
        {dataset_name}_mRNA_deg_{group_value}_intersect.csv
    """
    dataset_name = (
        validate_annotation_dataset_name(
            dataset_name
        )
    )
    group_value = (
        validate_scst_group_value_path_component(
            group_value
        )
    )

    return get_scst_group_value_result_path(
        annotation_dir=annotation_dir,
        filename=(
            f"{dataset_name}_mRNA_deg_"
            f"{group_value}_intersect.csv"
        ),
    )


def get_scst_dataset_log2fc_background_file_path(
    *,
    annotation_dir: Path,
    dataset_name: str,
    group_value: str,
) -> Path:
    dataset_name = (
        validate_annotation_dataset_name(
            dataset_name
        )
    )
    group_value = (
        validate_scst_group_value_path_component(
            group_value
        )
    )

    return get_scst_group_value_result_path(
        annotation_dir=annotation_dir,
        filename=(
            f"{dataset_name}_ceRNA_"
            f"{group_value}_background.csv"
        ),
    )


def get_scst_dataset_exp_correlation_file_path(
    *,
    annotation_dir: Path,
    dataset_name: str,
    group_value: str,
) -> Path:
    dataset_name = (
        validate_annotation_dataset_name(
            dataset_name
        )
    )
    group_value = (
        validate_scst_group_value_path_component(
            group_value
        )
    )

    return get_scst_group_value_result_path(
        annotation_dir=annotation_dir,
        filename=(
            f"{dataset_name}_ceRNA_corr_"
            f"{group_value}.csv"
        ),
    )


def get_scst_dataset_survival_file_path(
    *,
    annotation_dir: Path,
    dataset_name: str,
    group_value: str,
) -> Path:
    dataset_name = (
        validate_annotation_dataset_name(
            dataset_name
        )
    )
    group_value = (
        validate_scst_group_value_path_component(
            group_value
        )
    )

    return get_scst_group_value_result_path(
        annotation_dir=annotation_dir,
        filename=(
            f"{dataset_name}_survival_analysis_"
            f"{group_value}.csv"
        ),
    )


def get_scst_dataset_mrna_gsea_file_path(
    *,
    annotation_dir: Path,
    dataset_name: str,
    group_value: str,
) -> Path:
    dataset_name = (
        validate_annotation_dataset_name(
            dataset_name
        )
    )
    group_value = (
        validate_scst_group_value_path_component(
            group_value
        )
    )

    return get_scst_group_value_result_path(
        annotation_dir=annotation_dir,
        filename=(
            f"{dataset_name}_mRNA_gsea_"
            f"{group_value}.csv"
        ),
    )


def get_scst_dataset_cmdrug_result_dir_path(
    *,
    annotation_dir: Path,
    dataset_name: str,
    group_value: str,
) -> Path:
    dataset_name = (
        validate_annotation_dataset_name(
            dataset_name
        )
    )
    group_value = (
        validate_scst_group_value_path_component(
            group_value
        )
    )

    return get_scst_group_value_result_path(
        annotation_dir=annotation_dir,
        filename=(
            f"{dataset_name}_CMdrug_result_"
            f"{group_value}"
        ),
    )


def is_scst_existing_file(
    file_path: Path,
) -> bool:
    return (
        file_path.exists()
        and file_path.is_file()
    )


def is_scst_non_empty_file_directory(
    dir_path: Path,
) -> bool:
    """
    CMdrug is available only when the result directory exists
    and contains at least one direct child file.

    Nested subdirectories are intentionally not searched.
    """
    if (
        not dir_path.exists()
        or not dir_path.is_dir()
    ):
        return False

    try:
        return any(
            child.is_file()
            for child in dir_path.iterdir()
        )
    except OSError:
        return False


def build_scst_group_value_result_paths(
    *,
    annotation_dir: Path,
    dataset_name: str,
    group_value: str,
) -> dict:
    """
    Resolve all visualization resources for one SC/ST
    group value.
    """
    return {
        "cerna_axis_file": (
            get_scst_dataset_cerna_axis_file_path(
                annotation_dir=annotation_dir,
                dataset_name=dataset_name,
                group_value=group_value,
            )
        ),
        "map_immune_axis_file": (
            get_scst_dataset_map_immune_axis_file_path(
                annotation_dir=annotation_dir,
                dataset_name=dataset_name,
                group_value=group_value,
            )
        ),
        "axis_final_file": (
            get_scst_dataset_axis_final_file_path(
                annotation_dir=annotation_dir,
                dataset_name=dataset_name,
                group_value=group_value,
            )
        ),
        "cmap_file": (
            get_scst_dataset_cmap_file_path(
                annotation_dir=annotation_dir,
                dataset_name=dataset_name,
                group_value=group_value,
            )
        ),
        "deg_file": (
            get_scst_dataset_deg_file_path(
                annotation_dir=annotation_dir,
                dataset_name=dataset_name,
                group_value=group_value,
            )
        ),
        "deg_intersect_file": (
            get_scst_dataset_deg_intersect_file_path(
                annotation_dir=annotation_dir,
                dataset_name=dataset_name,
                group_value=group_value,
            )
        ),
        "log2fc_background_file": (
            get_scst_dataset_log2fc_background_file_path(
                annotation_dir=annotation_dir,
                dataset_name=dataset_name,
                group_value=group_value,
            )
        ),
        "exp_correlation_file": (
            get_scst_dataset_exp_correlation_file_path(
                annotation_dir=annotation_dir,
                dataset_name=dataset_name,
                group_value=group_value,
            )
        ),
        "survival_file": (
            get_scst_dataset_survival_file_path(
                annotation_dir=annotation_dir,
                dataset_name=dataset_name,
                group_value=group_value,
            )
        ),
        "mrna_gsea_file": (
            get_scst_dataset_mrna_gsea_file_path(
                annotation_dir=annotation_dir,
                dataset_name=dataset_name,
                group_value=group_value,
            )
        ),
        "cmdrug_result_dir": (
            get_scst_dataset_cmdrug_result_dir_path(
                annotation_dir=annotation_dir,
                dataset_name=dataset_name,
                group_value=group_value,
            )
        ),
    }
