from pathlib import Path

import pandas as pd


WORKFLOW_SPONGE_FILENAME_SUFFIX = "_sponge_result.csv"


# Sponge 文件中所有可能返回的字段。
WORKFLOW_SPONGE_COLUMNS = [
    "axis_id",
    "axis_type",
    "mRNA",
    "miRNA",
    "lncRNA",
    "circRNA",
    "cor",
    "pcor",
    "mscor",
]


# 所有 Sponge 结果都必须包含的核心字段。
WORKFLOW_SPONGE_REQUIRED_COLUMNS = [
    "axis_id",
    "axis_type",
    "mRNA",
    "miRNA",
    "cor",
    "pcor",
    "mscor",
]


WORKFLOW_SPONGE_NUMERIC_COLUMNS = {
    "cor",
    "pcor",
    "mscor",
}


class WorkflowSpongeInputError(ValueError):
    """Raised when a Sponge result file is invalid."""


class WorkflowSpongePathError(ValueError):
    """Raised when a Sponge result path is invalid."""


def validate_safe_name(value: str, field_name: str) -> None:
    """
    Validate a value that will be used as part of a filesystem path.
    """
    value = str(value or "").strip()

    if not value:
        raise WorkflowSpongePathError(
            f"Missing required parameter: {field_name}."
        )

    if "/" in value or "\\" in value or ".." in value:
        raise WorkflowSpongePathError(
            f"Invalid {field_name} parameter."
        )


def get_workflow_task_output_dir(task) -> Path:
    """
    Return the resolved output directory for a workflow task.
    """
    return Path(
        task.get_output_dir_absolute_path()
    ).resolve()


def get_workflow_sponge_file_path(task) -> Path:
    """
    Build the Sponge result file path.

    Expected filename:
        {task_name}_sponge_result.csv
    """
    task_name = str(task.task_name or "").strip()

    validate_safe_name(
        task_name,
        "task_name",
    )

    output_dir = get_workflow_task_output_dir(task)

    file_path = (
        output_dir
        / f"{task_name}{WORKFLOW_SPONGE_FILENAME_SUFFIX}"
    ).resolve()

    try:
        file_path.relative_to(output_dir)
    except ValueError as exc:
        raise WorkflowSpongePathError(
            "Invalid workflow Sponge result file path."
        ) from exc

    return file_path


def validate_workflow_sponge_file(task) -> Path:
    """
    Validate that the task Sponge result file exists.
    """
    file_path = get_workflow_sponge_file_path(task)

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            "Workflow Sponge result file not found: "
            f"{file_path.name}"
        )

    return file_path


def read_workflow_sponge_file(
    task,
    required_columns: list[str] | set[str] | None = None,
) -> tuple[Path, pd.DataFrame]:
    """
    Read and validate a task Sponge result file.
    """
    file_path = validate_workflow_sponge_file(task)

    return read_sponge_file_by_path(
        file_path=file_path,
        required_columns=(
            required_columns
            or WORKFLOW_SPONGE_REQUIRED_COLUMNS
        ),
    )


def validate_sponge_file_path(file_path: Path) -> Path:
    """
    Validate an explicitly provided Sponge file path.
    """
    file_path = Path(file_path).resolve()

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Sponge result file not found: {file_path.name}"
        )

    return file_path


def read_sponge_file_by_path(
    file_path: Path,
    required_columns: list[str] | set[str],
) -> tuple[Path, pd.DataFrame]:
    """
    Read a Sponge CSV and validate its required columns.
    """
    file_path = validate_sponge_file_path(file_path)

    try:
        df = pd.read_csv(file_path)
    except Exception as exc:
        raise WorkflowSpongeInputError(
            "Failed to read Sponge result file: "
            f"{str(exc)}"
        ) from exc

    required_columns = set(required_columns)
    existing_columns = set(df.columns)

    missing_columns = required_columns - existing_columns

    if missing_columns:
        raise WorkflowSpongeInputError(
            "Sponge result file is missing required column(s): "
            f"{', '.join(sorted(missing_columns))}."
        )

    return file_path, df


def resolve_existing_sponge_columns(
    *,
    df: pd.DataFrame,
    candidate_columns: list[str],
    required_columns: list[str] | set[str],
) -> list[str]:
    """
    Resolve ordered response columns.

    Rules:
    - All required columns must exist.
    - Optional columns are included only if present in the CSV.
    - Response order follows candidate_columns.
    """
    required_columns = set(required_columns)
    existing_columns = set(df.columns)

    missing_required_columns = (
        required_columns - existing_columns
    )

    if missing_required_columns:
        raise WorkflowSpongeInputError(
            "Sponge result file is missing required column(s): "
            f"{', '.join(sorted(missing_required_columns))}."
        )

    return [
        column
        for column in candidate_columns
        if column in existing_columns
    ]


def normalize_workflow_sponge_dataframe(
    df: pd.DataFrame,
    columns: list[str],
    numeric_columns: set[str] | None = None,
) -> pd.DataFrame:
    """
    Normalize Sponge dataframe columns and data types.
    """
    numeric_columns = (
        numeric_columns
        or WORKFLOW_SPONGE_NUMERIC_COLUMNS
    )

    normalized_df = df.copy()

    for column in columns:
        if column not in normalized_df.columns:
            normalized_df[column] = pd.NA

        if column in numeric_columns:
            normalized_df[column] = pd.to_numeric(
                normalized_df[column],
                errors="coerce",
            )
        else:
            normalized_df[column] = (
                normalized_df[column].astype("object")
            )

    normalized_df = normalized_df.replace(
        [float("inf"), float("-inf")],
        pd.NA,
    )

    return normalized_df[columns]


def safe_float_or_none(value):
    """
    Convert a value to a JSON-safe float.
    """
    if value is None or pd.isna(value):
        return None

    try:
        result = float(value)
    except (TypeError, ValueError):
        return None

    if pd.isna(result):
        return None

    return result


def safe_str_or_empty(value):
    """
    Convert a value to a JSON-safe string.
    """
    if value is None or pd.isna(value):
        return ""

    return str(value)


def serialize_workflow_sponge_dataframe(
    df: pd.DataFrame,
    columns: list[str],
    numeric_columns: set[str] | None = None,
) -> list[dict]:
    """
    Serialize a normalized Sponge dataframe into JSON-safe records.
    """
    numeric_columns = (
        numeric_columns
        or WORKFLOW_SPONGE_NUMERIC_COLUMNS
    )

    clean_df = df.replace(
        [float("inf"), float("-inf")],
        pd.NA,
    )

    records = []

    for _, row in clean_df.iterrows():
        item = {}

        for column in columns:
            value = row.get(column)

            if column in numeric_columns:
                item[column] = safe_float_or_none(value)
            else:
                item[column] = safe_str_or_empty(value)

        records.append(item)

    return records


def build_sponge_summary(df: pd.DataFrame) -> dict:
    """
    Build basic statistics for the Sponge result.

    This summary is optional but useful for the frontend control panel
    and result overview.
    """
    summary = {
        "raw_count": int(df.shape[0]),
        "axis_type_counts": {},
        "lncrna_axis_count": 0,
        "circrna_axis_count": 0,
    }

    if "axis_type" in df.columns:
        axis_type_series = (
            df["axis_type"]
            .dropna()
            .astype(str)
            .str.strip()
        )

        axis_type_series = axis_type_series[
            axis_type_series != ""
        ]

        summary["axis_type_counts"] = {
            str(axis_type): int(count)
            for axis_type, count
            in axis_type_series.value_counts().items()
        }

    if "lncRNA" in df.columns:
        lncrna_series = (
            df["lncRNA"]
            .dropna()
            .astype(str)
            .str.strip()
        )

        summary["lncrna_axis_count"] = int(
            (lncrna_series != "").sum()
        )

    if "circRNA" in df.columns:
        circrna_series = (
            df["circRNA"]
            .dropna()
            .astype(str)
            .str.strip()
        )

        summary["circrna_axis_count"] = int(
            (circrna_series != "").sum()
        )

    return summary


def build_sponge_response_from_dataframe(
    *,
    df: pd.DataFrame,
    sponge_file_name: str,
    columns: list[str] | None = None,
    required_columns: list[str] | set[str] | None = None,
    numeric_columns: set[str] | None = None,
    base_response: dict | None = None,
) -> dict:
    """
    Normalize and serialize a Sponge dataframe into an API response.
    """
    columns = columns or WORKFLOW_SPONGE_COLUMNS

    required_columns = (
        required_columns
        or WORKFLOW_SPONGE_REQUIRED_COLUMNS
    )

    numeric_columns = (
        numeric_columns
        or WORKFLOW_SPONGE_NUMERIC_COLUMNS
    )

    visible_columns = resolve_existing_sponge_columns(
        df=df,
        candidate_columns=columns,
        required_columns=required_columns,
    )

    normalized_df = normalize_workflow_sponge_dataframe(
        df=df,
        columns=visible_columns,
        numeric_columns=numeric_columns,
    )

    results = serialize_workflow_sponge_dataframe(
        df=normalized_df,
        columns=visible_columns,
        numeric_columns=numeric_columns,
    )

    response_data = {
        "sponge_file": sponge_file_name,
        "count": int(normalized_df.shape[0]),
        "columns": visible_columns,
        "summary": build_sponge_summary(
            normalized_df
        ),
        "results": results,
    }

    if base_response:
        response_data = {
            **base_response,
            **response_data,
        }

    return response_data
