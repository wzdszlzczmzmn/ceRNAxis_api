import math
from pathlib import Path

import numpy as np
import pandas as pd


CM_SCORE_RESULT_DIR_NAME = "CM_results"
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


def get_workflow_cm_results_dir(task) -> Path:
    """
    Return:
        task_output_dir / CM_results
    """
    output_dir = get_workflow_task_output_dir(task)

    cm_results_dir = (
        output_dir / CM_SCORE_RESULT_DIR_NAME
    ).resolve()

    try:
        cm_results_dir.relative_to(output_dir)
    except ValueError as exc:
        raise WorkflowCMScorePathError(
            "Invalid CM-results directory path."
        ) from exc

    return cm_results_dir


def validate_workflow_cm_results_dir(task) -> Path:
    cm_results_dir = get_workflow_cm_results_dir(task)

    if not cm_results_dir.exists():
        raise FileNotFoundError(
            f"CM-results directory not found: "
            f"{CM_SCORE_RESULT_DIR_NAME}"
        )

    if not cm_results_dir.is_dir():
        raise WorkflowCMScorePathError(
            f"CM-results path is not a directory: "
            f"{CM_SCORE_RESULT_DIR_NAME}"
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


def get_workflow_cm_score_file_path(
    task,
    item_value: str,
) -> Path:
    cm_results_dir = get_workflow_cm_results_dir(
        task
    )

    filename = get_workflow_cm_score_filename(
        item_value
    )

    file_path = (
        cm_results_dir / filename
    ).resolve()

    try:
        file_path.relative_to(cm_results_dir)
    except ValueError as exc:
        raise WorkflowCMScorePathError(
            "Invalid CM-score result file path."
        ) from exc

    return file_path


def validate_workflow_cm_score_file(
    task,
    item_value: str,
) -> Path:
    file_path = get_workflow_cm_score_file_path(
        task=task,
        item_value=item_value,
    )

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"CM-score result file not found: "
            f"{file_path.name}"
        )

    return file_path


def get_available_workflow_cm_score_items(
    task,
) -> list[dict]:
    """
    Return all valid CM-score files directly inside:

        task_output_dir/CM_results/

    Only files matching:
        {item}_CM_scores.csv

    are returned.

    The filesystem path is deliberately not exposed to the frontend.
    """
    cm_results_dir = get_workflow_cm_results_dir(
        task
    )

    if (
        not cm_results_dir.exists()
        or not cm_results_dir.is_dir()
    ):
        return []

    items = []

    for file_path in cm_results_dir.iterdir():
        if not file_path.is_file():
            continue

        item_value = extract_cm_score_item_value(
            file_path.name
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

    items.sort(
        key=lambda item: (
            item["label"].casefold(),
            item["label"],
        )
    )

    return items


def build_workflow_cm_score_options_response(
    *,
    task,
    task_type: str,
) -> dict:
    items = get_available_workflow_cm_score_items(
        task
    )

    return {
        "uuid": str(task.uuid),
        "task_type": task_type,
        "task_name": task.task_name,
        "count": len(items),
        "default_item": (
            items[0]["value"]
            if items
            else None
        ),
        "results": items,
    }


def read_workflow_cm_score_file(
    task,
    item_value: str,
) -> tuple[Path, pd.DataFrame]:
    file_path = validate_workflow_cm_score_file(
        task=task,
        item_value=item_value,
    )

    try:
        dataframe = pd.read_csv(file_path)
    except UnicodeDecodeError as exc:
        raise WorkflowCMScoreInputError(
            f"CM-score file must be UTF-8 encoded: {file_path.name}."
        ) from exc
    except Exception as exc:
        raise WorkflowCMScoreInputError(
            f"Failed to read CM-score result file: {exc}"
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


def build_workflow_cm_score_result_response(
    *,
    task,
    task_type: str,
    item_value: str,
    file_path: Path,
    dataframe: pd.DataFrame,
) -> dict:
    results = serialize_workflow_cm_score_dataframe(
        dataframe
    )

    dataset_options = get_cm_score_dataset_options(
        results
    )

    return {
        "uuid": str(task.uuid),
        "task_type": task_type,
        "task_name": task.task_name,

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
            "intercept_key": "a",
            "intercept": 0.4986,
            "score_field": "cm_score",
            "score_label": "CM-Score",
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
