import math
from pathlib import Path

import numpy as np
import pandas as pd


CUSTOM_LIST_CM_SCORE_RESULT_DIR_NAME = "CM_results"
CMDRUG_RESULT_DIR_SUFFIX = "_CMdrug_result"
CM_SCORE_FILE_SUFFIX = "_CM_scores.csv"


CM_SCORE_DATASET_COLUMN = "dataset"
CM_SCORE_SCORE_COLUMN = "CM_Score"

CM_SCORE_PATHWAY_FIELDS = [
    "Antigen_Processing_and_Presentation",
    "NaturalKiller_Cell_Cytotoxicity",
    "TCR_Signaling_Pathway",
    "Cytotoxiclty_of_ImmuCellAI",
    "Antimicrobials",
    "BCR_Signaling_Pathway",
]

CM_SCORE_REQUIRED_COLUMNS = {
    CM_SCORE_DATASET_COLUMN,
    CM_SCORE_SCORE_COLUMN,
    *CM_SCORE_PATHWAY_FIELDS,
}

CM_SCORE_OPTIONAL_COLUMNS = [
    "pass_relaxed_all_padj_lt_0.2_NES_gt_0",
    "strict_pathway_count_padj_lt_0.05_NES_gt_0",
    "pass_strategy2",
    "cell_line",
]

CM_SCORE_FIELD_CONFIG = [
    {
        "source_field": "Antigen_Processing_and_Presentation",
        "response_field": "antigen_processing_and_presentation",
        "label": "CMDrug C1: Antigen Processing and Presentation",
        "short_label": "C1",
        "weight_key": "b",
        "weight": 0.0969,
    },
    {
        "source_field": "NaturalKiller_Cell_Cytotoxicity",
        "response_field": "natural_killer_cell_cytotoxicity",
        "label": "CMDrug C2: Natural Killer Cell Cytotoxicity",
        "short_label": "C2",
        "weight_key": "c",
        "weight": 0.0969,
    },
    {
        "source_field": "TCR_Signaling_Pathway",
        "response_field": "tcr_signaling_pathway",
        "label": "CMDrug C3: TCR Signaling Pathway",
        "short_label": "C3",
        "weight_key": "d",
        "weight": 0.0307,
    },
    {
        "source_field": "Cytotoxiclty_of_ImmuCellAI",
        "response_field": "cytotoxicity_of_immucellai",
        "label": "CMDrug C4: Cytotoxicity of ImmuCellAI",
        "short_label": "C4",
        "weight_key": "e",
        "weight": 0.0117,
    },
    {
        "source_field": "Antimicrobials",
        "response_field": "antimicrobials",
        "label": "CMDrug M1: Antimicrobials",
        "short_label": "M1",
        "weight_key": "f",
        "weight": 0.0124,
    },
    {
        "source_field": "BCR_Signaling_Pathway",
        "response_field": "bcr_signaling_pathway",
        "label": "CMDrug M2: BCR Signaling Pathway",
        "short_label": "M2",
        "weight_key": "g",
        "weight": 0.0213,
    },
]

CM_SCORE_FORMULA_CONFIG = {
    "intercept_key": "a",
    "intercept": 0.4986,
    "score_field": "CM_Score",
    "score_label": "CM-Score",
}


class WorkflowCMScoreInputError(ValueError):
    pass


class WorkflowCMScorePathError(ValueError):
    pass


def get_workflow_task_output_dir(task) -> Path:
    """
    Return the resolved workflow task output directory.

    CustomListQueryTask, PairedCohortTask and HybridReferenceTask
    all expose get_output_dir_absolute_path().
    """
    if not hasattr(task, "get_output_dir_absolute_path"):
        raise WorkflowCMScorePathError(
            "Task does not provide an output directory."
        )

    output_dir = Path(
        task.get_output_dir_absolute_path()
    ).resolve()

    return output_dir


def get_workflow_cm_results_dir(
    task,
    *,
    group_value: str | None = None,
) -> Path:
    output_dir = get_workflow_task_output_dir(task)

    task_type = task.__class__.__name__
    task_name = str(
        getattr(task, "task_name", "") or ""
    ).strip()

    if not task_name:
        raise WorkflowCMScorePathError(
            "Task does not define task_name."
        )

    if (
        "/" in task_name
        or "\\" in task_name
        or ".." in task_name
    ):
        raise WorkflowCMScorePathError(
            "Invalid task_name."
        )

    if task_type == "CustomListQueryTask":
        dir_name = CUSTOM_LIST_CM_SCORE_RESULT_DIR_NAME

    elif task_type in {
        "PairedCohortTask",
        "HybridReferenceTask",
    }:
        dir_name = (
            f"{task_name}"
            f"{CMDRUG_RESULT_DIR_SUFFIX}"
        )

    elif task_type == "SCSTHybridReferenceTask":
        group_value = str(
            group_value or ""
        ).strip()

        if not group_value:
            raise WorkflowCMScoreInputError(
                "Missing required parameter: groupValue."
            )

        if (
            "/" in group_value
            or "\\" in group_value
            or ".." in group_value
        ):
            raise WorkflowCMScoreInputError(
                "Invalid groupValue parameter."
            )

        dir_name = (
            f"{task_name}"
            f"{CMDRUG_RESULT_DIR_SUFFIX}_"
            f"{group_value}"
        )

    else:
        raise WorkflowCMScorePathError(
            f"Unsupported CM-score task type: {task_type}."
        )

    cm_results_dir = (
        output_dir / dir_name
    ).resolve()

    try:
        cm_results_dir.relative_to(output_dir)
    except ValueError as exc:
        raise WorkflowCMScorePathError(
            "Invalid CM-results directory path."
        ) from exc

    return cm_results_dir


def validate_workflow_cm_results_dir(
    task,
    *,
    group_value: str | None = None,
) -> Path:
    cm_results_dir = get_workflow_cm_results_dir(
        task,
        group_value=group_value,
    )

    if not cm_results_dir.exists():
        raise FileNotFoundError(
            "CM-results directory not found: "
            f"{cm_results_dir.name}"
        )

    if not cm_results_dir.is_dir():
        raise WorkflowCMScorePathError(
            "CM-results path is not a directory: "
            f"{cm_results_dir.name}"
        )

    return cm_results_dir


def extract_cm_score_item_value(
    filename: str,
) -> str | None:
    """
    Parse:
        {item}_CM_scores.csv

    Examples:
        TP53_CM_scores.csv
            -> TP53

        axis_000001_CM_scores.csv
            -> axis_000001
    """
    filename = str(filename or "").strip()

    if not filename.endswith(CM_SCORE_FILE_SUFFIX):
        return None

    item_value = filename[
        :-len(CM_SCORE_FILE_SUFFIX)
    ].strip()

    if not item_value:
        return None

    return item_value


def validate_cm_score_item_value(
    item_value: str,
) -> str:
    """
    Validate a value that will later be used to locate a CM-score file.

    The value may be a gene name or an axis ID, so validation should not
    impose biological semantics. It only prevents path traversal.
    """
    normalized_value = str(
        item_value or ""
    ).strip()

    if not normalized_value:
        raise WorkflowCMScoreInputError(
            "Missing required parameter: item."
        )

    if (
        "/" in normalized_value
        or "\\" in normalized_value
        or ".." in normalized_value
    ):
        raise WorkflowCMScoreInputError(
            "Invalid item parameter."
        )

    if normalized_value in {".", ".."}:
        raise WorkflowCMScoreInputError(
            "Invalid item parameter."
        )

    return normalized_value


def get_workflow_cm_score_filename(
    item_value: str,
) -> str:
    normalized_value = validate_cm_score_item_value(
        item_value
    )

    return (
        f"{normalized_value}"
        f"{CM_SCORE_FILE_SUFFIX}"
    )


def get_cm_score_file_path_from_dir(
    *,
    cm_results_dir: Path,
    item_value: str,
) -> Path:
    """
    Resolve one {item}_CM_scores.csv path inside an already-resolved
    CM-results directory.

    This helper is intentionally task/source independent so Workflow
    and Dataset Annotation can share the same item/path validation.
    """
    cm_results_dir = Path(
        cm_results_dir
    ).resolve()

    filename = get_workflow_cm_score_filename(
        item_value
    )

    file_path = (
        cm_results_dir
        / filename
    ).resolve()

    try:
        file_path.relative_to(
            cm_results_dir
        )
    except ValueError as exc:
        raise WorkflowCMScorePathError(
            "Invalid CM-score result file path."
        ) from exc

    return file_path


def get_workflow_cm_score_file_path(
    task,
    item_value: str,
    *,
    group_value: str | None = None,
) -> Path:
    cm_results_dir = get_workflow_cm_results_dir(
        task,
        group_value=group_value,
    )

    return get_cm_score_file_path_from_dir(
        cm_results_dir=cm_results_dir,
        item_value=item_value,
    )


def validate_workflow_cm_score_file(
    task,
    item_value: str,
    *,
    group_value: str | None = None,
) -> Path:
    file_path = get_workflow_cm_score_file_path(
        task=task,
        item_value=item_value,
        group_value=group_value,
    )

    if (
        not file_path.exists()
        or not file_path.is_file()
    ):
        raise FileNotFoundError(
            "CM-score result file not found: "
            f"{file_path.name}"
        )

    return file_path


def get_available_cm_score_items_from_dir(
    *,
    cm_results_dir: Path,
) -> list[dict]:
    """
    Scan one CM-results directory for {item}_CM_scores.csv files.

    Missing/non-directory paths intentionally return an empty option
    list. This preserves the existing Workflow Options endpoint
    behavior and is also suitable for Dataset Annotation.
    """
    cm_results_dir = Path(
        cm_results_dir
    ).resolve()

    if (
        not cm_results_dir.exists()
        or not cm_results_dir.is_dir()
    ):
        return []

    items = []

    try:
        children = cm_results_dir.iterdir()

        for file_path in children:
            if not file_path.is_file():
                continue

            item_value = (
                extract_cm_score_item_value(
                    file_path.name
                )
            )

            if item_value is None:
                continue

            items.append(
                {
                    "value": item_value,
                    "label": item_value,
                    "file_name": file_path.name,
                }
            )

    except OSError as exc:
        raise WorkflowCMScorePathError(
            "Failed to inspect CM-results directory."
        ) from exc

    items.sort(
        key=lambda item: (
            item["label"].casefold(),
            item["label"],
        )
    )

    return items


def get_available_workflow_cm_score_items(
    task,
    *,
    group_value: str | None = None,
) -> list[dict]:
    cm_results_dir = (
        get_workflow_cm_results_dir(
            task,
            group_value=group_value,
        )
    )

    return get_available_cm_score_items_from_dir(
        cm_results_dir=cm_results_dir,
    )


def build_cm_score_options_response_common(
    *,
    base_response: dict | None,
    items: list[dict],
) -> dict:
    response = {
        "count": len(items),
        "default_item": (
            items[0]["value"]
            if items
            else None
        ),
        "results": items,
    }

    if base_response:
        response = {
            **base_response,
            **response,
        }

    return response


def build_workflow_cm_score_options_response(
    *,
    task,
    task_type: str,
    group_value: str | None = None,
) -> dict:
    items = get_available_workflow_cm_score_items(
        task,
        group_value=group_value,
    )

    base_response = {
        "uuid": str(task.uuid),
        "task_type": task_type,
        "task_name": task.task_name,
    }

    if group_value is not None:
        base_response["group_value"] = group_value

    return build_cm_score_options_response_common(
        base_response=base_response,
        items=items,
    )


def read_cm_score_file_by_path(
    *,
    file_path: Path,
) -> tuple[Path, pd.DataFrame]:
    """
    Read and validate a CM-score CSV from an already-resolved path.
    """
    file_path = Path(
        file_path
    ).resolve()

    if (
        not file_path.exists()
        or not file_path.is_file()
    ):
        raise FileNotFoundError(
            "CM-score result file not found: "
            f"{file_path.name}"
        )

    try:
        dataframe = pd.read_csv(
            file_path
        )
    except UnicodeDecodeError as exc:
        raise WorkflowCMScoreInputError(
            "CM-score file must be UTF-8 encoded: "
            f"{file_path.name}."
        ) from exc
    except Exception as exc:
        raise WorkflowCMScoreInputError(
            "Failed to read CM-score result file: "
            f"{exc}"
        ) from exc

    missing_columns = (
        CM_SCORE_REQUIRED_COLUMNS
        - set(dataframe.columns)
    )

    if missing_columns:
        raise WorkflowCMScoreInputError(
            "CM-score result file is missing required columns: "
            f"{', '.join(sorted(missing_columns))}."
        )

    return file_path, dataframe


def read_workflow_cm_score_file(
    task,
    item_value: str,
    *,
    group_value: str | None = None,
) -> tuple[Path, pd.DataFrame]:
    file_path = (
        validate_workflow_cm_score_file(
            task=task,
            item_value=item_value,
            group_value=group_value,
        )
    )

    return read_cm_score_file_by_path(
        file_path=file_path,
    )


def normalize_optional_float(value):
    if value is None:
        return None

    if isinstance(value, np.generic):
        value = value.item()

    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return None

    if not math.isfinite(numeric_value):
        return None

    return numeric_value


def normalize_optional_bool(value):
    if value is None:
        return None

    if isinstance(value, np.bool_):
        return bool(value)

    if isinstance(value, bool):
        return value

    normalized = str(value).strip().lower()

    if normalized in {"true", "1", "yes"}:
        return True

    if normalized in {"false", "0", "no"}:
        return False

    return None


def normalize_optional_int(value):
    if value is None:
        return None

    try:
        numeric_value = int(value)
    except (TypeError, ValueError):
        return None

    return numeric_value


def serialize_workflow_cm_score_dataframe(
    dataframe: pd.DataFrame,
) -> list[dict]:
    results = []

    for row_number, (_, row) in enumerate(
        dataframe.iterrows(),
        start=1,
    ):
        dataset = str(
            row.get("dataset", "") or ""
        ).strip()

        if not dataset:
            continue

        item = {
            "id": row_number,
            "dataset": dataset,
            "cm_score": normalize_optional_float(
                row.get("CM_Score")
            ),
            "cell_line": str(
                row.get("cell_line", "") or ""
            ).strip() or None,
            "pass_relaxed_all_padj_lt_0_2_nes_gt_0": (
                normalize_optional_bool(
                    row.get(
                        "pass_relaxed_all_padj_lt_0.2_NES_gt_0"
                    )
                )
            ),
            "strict_pathway_count_padj_lt_0_05_nes_gt_0": (
                normalize_optional_int(
                    row.get(
                        "strict_pathway_count_padj_lt_0.05_NES_gt_0"
                    )
                )
            ),
            "pass_strategy2": normalize_optional_bool(
                row.get("pass_strategy2")
            ),
            "pathway_values": {},
        }

        for config in CM_SCORE_FIELD_CONFIG:
            item["pathway_values"][
                config["response_field"]
            ] = normalize_optional_float(
                row.get(config["source_field"])
            )

        results.append(item)

    return results


def get_cm_score_dataset_options(
    results: list[dict],
) -> list[dict]:
    seen = set()
    options = []

    for row in results:
        dataset = str(
            row.get("dataset", "")
            or ""
        ).strip()

        if not dataset or dataset in seen:
            continue

        seen.add(dataset)

        options.append(
            {
                "value": dataset,
                "label": dataset,
            }
        )

    return options


def build_cm_score_result_response_common(
    *,
    base_response: dict | None,
    item_value: str,
    file_path: Path,
    dataframe: pd.DataFrame,
) -> dict:
    """
    Build the source-independent CM-score response consumed by the
    shared frontend visualization.
    """
    item_value = validate_cm_score_item_value(
        item_value
    )

    results = serialize_workflow_cm_score_dataframe(
        dataframe
    )

    dataset_options = get_cm_score_dataset_options(
        results
    )

    response = {
        "item": item_value,
        "cm_score_file": file_path.name,

        "count": len(results),
        "default_dataset": (
            dataset_options[0]["value"]
            if dataset_options
            else None
        ),
        "dataset_options": dataset_options,

        "plot": {
            "dataset_field": "dataset",
            "score_field": "cm_score",
            "pathway_values_field": "pathway_values",
            "positive_color": "red",
            "negative_color": "blue",
            "point_size_mode": "constant",
        },

        "formula": {
            "intercept_key": (
                CM_SCORE_FORMULA_CONFIG[
                    "intercept_key"
                ]
            ),
            "intercept": (
                CM_SCORE_FORMULA_CONFIG[
                    "intercept"
                ]
            ),
            "score_field": "cm_score",
            "score_label": (
                CM_SCORE_FORMULA_CONFIG[
                    "score_label"
                ]
            ),
            "components": [
                {
                    "field": config["response_field"],
                    "source_field": config["source_field"],
                    "label": config["label"],
                    "short_label": config["short_label"],
                    "weight_key": config["weight_key"],
                    "weight": config["weight"],
                }
                for config in CM_SCORE_FIELD_CONFIG
            ],
        },

        "results": results,
    }

    if base_response:
        response = {
            **base_response,
            **response,
        }

    return response


def build_workflow_cm_score_result_response(
    *,
    task,
    task_type: str,
    item_value: str,
    file_path: Path,
    dataframe: pd.DataFrame,
    group_value: str | None = None,
) -> dict:
    base_response = {
        "uuid": str(task.uuid),
        "task_type": task_type,
        "task_name": task.task_name,
    }

    if group_value is not None:
        base_response["group_value"] = group_value

    return build_cm_score_result_response_common(
        base_response=base_response,
        item_value=item_value,
        file_path=file_path,
        dataframe=dataframe,
    )
