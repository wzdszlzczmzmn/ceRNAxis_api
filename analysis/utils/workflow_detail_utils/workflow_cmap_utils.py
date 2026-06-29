from pathlib import Path

import numpy as np
import pandas as pd


WORKFLOW_CMAP_FILENAME_SUFFIX = "_CMap.csv"

WORKFLOW_CMAP_REQUIRED_COLUMNS = set()

HYBRID_REFERENCE_CMAP_COLUMNS = [
    "c_perturbation",
    "c_perturbation_name",
    "c_perturbation_type",
    "n_tau",
    "c_cell_line",
    "n_perturbation_size",
]

HYBRID_REFERENCE_CMAP_REQUIRED_COLUMNS = set(
    HYBRID_REFERENCE_CMAP_COLUMNS
)


class WorkflowCMapInputError(ValueError):
    pass


class WorkflowCMapPathError(ValueError):
    pass


def validate_safe_name(value: str, field_name: str) -> None:
    value = str(value or "").strip()

    if not value:
        raise WorkflowCMapPathError(
            f"Missing required parameter: {field_name}."
        )

    if "/" in value or "\\" in value or ".." in value:
        raise WorkflowCMapPathError(
            f"Invalid {field_name} parameter."
        )


def get_workflow_task_output_dir(task) -> Path:
    return Path(task.get_output_dir_absolute_path()).resolve()


def get_workflow_cmap_file_path(task) -> Path:
    task_name = str(task.task_name).strip()

    validate_safe_name(task_name, "task_name")

    output_dir = get_workflow_task_output_dir(task)

    file_path = (
        output_dir / f"{task_name}{WORKFLOW_CMAP_FILENAME_SUFFIX}"
    ).resolve()

    if not str(file_path).startswith(str(output_dir)):
        raise WorkflowCMapPathError(
            "Invalid workflow CMap result file path."
        )

    return file_path


def validate_workflow_cmap_file(task) -> Path:
    file_path = get_workflow_cmap_file_path(task)

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Workflow CMap result file not found: {file_path.name}"
        )

    return file_path


def normalize_workflow_cmap_dataframe(
    df: pd.DataFrame,
) -> pd.DataFrame:
    normalized_df = df.copy()

    normalized_df = normalized_df.replace(
        [np.inf, -np.inf],
        np.nan,
    )

    normalized_df = normalized_df.where(
        pd.notnull(normalized_df),
        None,
    )

    return normalized_df


def serialize_workflow_cmap_dataframe(
    df: pd.DataFrame,
) -> list[dict]:
    normalized_df = normalize_workflow_cmap_dataframe(df)

    records = normalized_df.to_dict(orient="records")

    # Convert numpy scalar values into Python native values.
    serialized_records = []

    for record in records:
        item = {}

        for key, value in record.items():
            if isinstance(value, np.generic):
                value = value.item()

            item[str(key)] = value

        serialized_records.append(item)

    return serialized_records


def validate_cmap_file_path(file_path: Path) -> Path:
    file_path = Path(file_path).resolve()

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"CMap result file not found: {file_path.name}"
        )

    return file_path


def read_cmap_file_by_path(
    file_path: Path,
    required_columns: set[str] | list[str] | None = None,
) -> tuple[Path, pd.DataFrame]:
    file_path = validate_cmap_file_path(file_path)

    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        raise WorkflowCMapInputError(
            f"Failed to read CMap result file: {str(e)}"
        )

    required_columns = set(required_columns or [])

    if required_columns:
        missing_columns = required_columns - set(df.columns)

        if missing_columns:
            raise WorkflowCMapInputError(
                "CMap result file is missing required column(s): "
                f"{', '.join(sorted(missing_columns))}."
            )

    return file_path, df


def read_workflow_cmap_file(
    task,
    required_columns: set[str] | list[str] | None = None,
) -> tuple[Path, pd.DataFrame]:
    file_path = validate_workflow_cmap_file(task)

    return read_cmap_file_by_path(
        file_path=file_path,
        required_columns=required_columns,
    )


def build_cmap_response_from_dataframe(
    *,
    df: pd.DataFrame,
    cmap_file_name: str,
    base_response: dict | None = None,
) -> dict:
    normalized_df = normalize_workflow_cmap_dataframe(df)

    columns = [
        str(col)
        for col in normalized_df.columns
    ]

    results = serialize_workflow_cmap_dataframe(normalized_df)

    response_data = {
        "cmap_file": cmap_file_name,
        "columns": columns,
        "summary": {
            "raw_count": len(results),
        },
        "results": results,
    }

    if base_response:
        response_data = {
            **base_response,
            **response_data,
        }

    return response_data
