from pathlib import Path

import numpy as np
import pandas as pd


WORKFLOW_LOG2FC_BACKGROUND_FILENAME_SUFFIX = "_ceRNA_background.csv"

WORKFLOW_LOG2FC_BACKGROUND_REQUIRED_COLUMNS = {
    "miRNA",
    "ceRNA",
    "species",
    "database",
    "type",
    "miRNA_log2FC",
    "ceRNA_log2FC",
}

WORKFLOW_LOG2FC_BACKGROUND_STRING_COLUMNS = [
    "miRNA",
    "ceRNA",
    "species",
    "database",
    "type",
]

WORKFLOW_LOG2FC_BACKGROUND_NUMERIC_COLUMNS = [
    "miRNA_log2FC",
    "ceRNA_log2FC",
]

PAIRED_COHORT_VALID_BACKGROUND_TYPES = [
    "miRNA-mRNA",
    "miRNA-lncRNA",
    "miRNA-circRNA",
]

HYBRID_REFERENCE_VALID_BACKGROUND_TYPES = [
    "miRNA-mRNA",
    "miRNA-lncRNA",
]

WORKFLOW_LOG2FC_X_COL = "ceRNA_log2FC"
WORKFLOW_LOG2FC_Y_COL = "miRNA_log2FC"


class WorkflowLog2FCBackgroundInputError(ValueError):
    pass


class WorkflowLog2FCBackgroundPathError(ValueError):
    pass


def validate_safe_name(value: str, field_name: str) -> None:
    value = str(value or "").strip()

    if not value:
        raise WorkflowLog2FCBackgroundPathError(
            f"Missing required parameter: {field_name}."
        )

    if "/" in value or "\\" in value or ".." in value:
        raise WorkflowLog2FCBackgroundPathError(
            f"Invalid {field_name} parameter."
        )


def get_workflow_task_output_dir(task) -> Path:
    return Path(task.get_output_dir_absolute_path()).resolve()


def get_workflow_log2fc_background_file_path(task) -> Path:
    task_name = str(task.task_name).strip()

    validate_safe_name(task_name, "task_name")

    output_dir = get_workflow_task_output_dir(task)

    file_path = (
        output_dir / f"{task_name}{WORKFLOW_LOG2FC_BACKGROUND_FILENAME_SUFFIX}"
    ).resolve()

    if not str(file_path).startswith(str(output_dir)):
        raise WorkflowLog2FCBackgroundPathError(
            "Invalid workflow ceRNA background file path."
        )

    return file_path


def validate_workflow_log2fc_background_file(task) -> Path:
    file_path = get_workflow_log2fc_background_file_path(task)

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Workflow ceRNA background file not found: {file_path.name}"
        )

    return file_path


def read_workflow_log2fc_background_file(task) -> tuple[Path, pd.DataFrame]:
    file_path = validate_workflow_log2fc_background_file(task)

    return read_log2fc_background_file_by_path(file_path)


def normalize_workflow_log2fc_background_dataframe(
    df: pd.DataFrame,
) -> pd.DataFrame:
    normalized_df = df.copy()

    for col in WORKFLOW_LOG2FC_BACKGROUND_STRING_COLUMNS:
        if col in normalized_df.columns:
            normalized_df[col] = (
                normalized_df[col]
                .astype(str)
                .str.strip()
            )

    for col in WORKFLOW_LOG2FC_BACKGROUND_NUMERIC_COLUMNS:
        if col in normalized_df.columns:
            normalized_df[col] = pd.to_numeric(
                normalized_df[col],
                errors="coerce",
            )

    normalized_df = normalized_df.replace(
        [np.inf, -np.inf],
        np.nan,
    )

    return normalized_df


def get_workflow_available_background_types(
    df: pd.DataFrame,
    valid_types: list[str],
) -> list[str]:
    normalized_df = normalize_workflow_log2fc_background_dataframe(df)

    observed_types = set(
        normalized_df["type"]
        .dropna()
        .astype(str)
        .str.strip()
        .tolist()
    )

    return [
        type_value
        for type_value in valid_types
        if type_value in observed_types
    ]


def serialize_log2fc_background_points(
    df: pd.DataFrame,
    x_col: str = WORKFLOW_LOG2FC_X_COL,
    y_col: str = WORKFLOW_LOG2FC_Y_COL,
) -> list[dict]:
    points = []

    for _, row in df.iterrows():
        species = row.get("species", "")
        database = row.get("database", "")

        if pd.isna(species):
            species = ""

        if pd.isna(database):
            database = ""

        points.append(
            {
                "miRNA": row["miRNA"],
                "ceRNA": row["ceRNA"],
                "species": species,
                "database": database,
                "type": row["type"],
                "ceRNA_log2FC": float(row[x_col]),
                "miRNA_log2FC": float(row[y_col]),
                "anti_correlation": bool(row["anti_correlation"]),
            }
        )

    return points


def build_workflow_log2fc_correlation_response_data(
    *,
    task,
    task_type: str,
    df: pd.DataFrame,
    type_value: str,
    background_file_name: str,
    available_types: list[str],
    x_col: str = WORKFLOW_LOG2FC_X_COL,
    y_col: str = WORKFLOW_LOG2FC_Y_COL,
) -> dict:
    base_response = {
        "uuid": str(task.uuid),
        "task_type": task_type,
        "task_name": task.task_name,
    }

    return build_log2fc_correlation_response_data_from_dataframe(
        df=df,
        type_value=type_value,
        background_file_name=background_file_name,
        available_types=available_types,
        base_response=base_response,
        x_col=x_col,
        y_col=y_col,
    )


def get_available_workflow_log2fc_background_types(
    task,
    valid_types: list[str],
) -> list[str]:
    try:
        _, df = read_workflow_log2fc_background_file(task)
    except (
        FileNotFoundError,
        WorkflowLog2FCBackgroundPathError,
        WorkflowLog2FCBackgroundInputError,
    ):
        return []

    return get_workflow_available_background_types(
        df=df,
        valid_types=valid_types,
    )


def validate_log2fc_background_file_path(file_path: Path) -> Path:
    file_path = Path(file_path).resolve()

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"ceRNA background file not found: {file_path.name}"
        )

    return file_path


def read_log2fc_background_file_by_path(
    file_path: Path,
) -> tuple[Path, pd.DataFrame]:
    file_path = validate_log2fc_background_file_path(file_path)

    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        raise WorkflowLog2FCBackgroundInputError(
            f"Failed to read ceRNA background file: {str(e)}"
        )

    missing_columns = (
        WORKFLOW_LOG2FC_BACKGROUND_REQUIRED_COLUMNS
        - set(df.columns)
    )

    if missing_columns:
        raise WorkflowLog2FCBackgroundInputError(
            "ceRNA background file is missing required column(s): "
            f"{', '.join(sorted(missing_columns))}."
        )

    return file_path, df


def build_log2fc_correlation_response_data_from_dataframe(
    *,
    df: pd.DataFrame,
    type_value: str,
    background_file_name: str,
    available_types: list[str],
    base_response: dict | None = None,
    x_col: str = WORKFLOW_LOG2FC_X_COL,
    y_col: str = WORKFLOW_LOG2FC_Y_COL,
) -> dict:
    df = normalize_workflow_log2fc_background_dataframe(df)

    df = df[df["type"] == type_value].copy()

    raw_count = int(df.shape[0])

    required_subset = [
        "miRNA",
        "ceRNA",
        "type",
        x_col,
        y_col,
    ]

    df = (
        df
        .replace([np.inf, -np.inf], np.nan)
        .dropna(subset=required_subset)
        .copy()
    )

    df = df[
        (df["miRNA"].astype(str).str.strip() != "")
        & (df["ceRNA"].astype(str).str.strip() != "")
        & (df["type"].astype(str).str.strip() != "")
    ].copy()

    cleaned_count = int(df.shape[0])
    dropped_count = raw_count - cleaned_count

    response_data = {
        "type": type_value,
        "available_types": available_types,
        "background_file": background_file_name,
        "summary": {
            "raw_count": raw_count,
            "cleaned_count": cleaned_count,
            "dropped_count": dropped_count,
            "anti_count": 0,
            "same_count": 0,
        },
        "points": [],
    }

    if cleaned_count > 0:
        df["anti_correlation"] = df[x_col] * df[y_col] < 0

        anti_count = int(df["anti_correlation"].sum())
        same_count = cleaned_count - anti_count

        response_data["summary"].update(
            {
                "anti_count": anti_count,
                "same_count": same_count,
            }
        )

        response_data["points"] = serialize_log2fc_background_points(
            df=df,
            x_col=x_col,
            y_col=y_col,
        )

    if base_response:
        response_data = {
            **base_response,
            **response_data,
        }

    return response_data
