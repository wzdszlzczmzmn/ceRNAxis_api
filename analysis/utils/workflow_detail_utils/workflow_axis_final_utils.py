from pathlib import Path

import pandas as pd


WORKFLOW_AXIS_FINAL_FILENAME_SUFFIX = "_ceRNA_axis_final.csv"

WORKFLOW_AXIS_FINAL_FILENAME_SUFFIX = "_ceRNA_axis_final.csv"

WORKFLOW_AXIS_FINAL_CORE_COLUMNS = [
    "axis_id",
    "axis_regulation",
    "axis_type",
    "mRNA",
    "mRNA_log2FC",
    "mRNA_regulation",
    "miRNA",
    "miRNA_log2FC",
    "miRNA_regulation",
]

WORKFLOW_AXIS_FINAL_LNCRNA_COLUMNS = [
    "lncRNA",
    "lncRNA_log2FC",
    "lncRNA_regulation",
]

WORKFLOW_AXIS_FINAL_CIRCRNA_COLUMNS = [
    "circRNA",
    "circRNA_log2FC",
    "circRNA_regulation",
]

PAIRED_COHORT_AXIS_FINAL_COLUMNS = [
    *WORKFLOW_AXIS_FINAL_CORE_COLUMNS,
    *WORKFLOW_AXIS_FINAL_LNCRNA_COLUMNS,
    *WORKFLOW_AXIS_FINAL_CIRCRNA_COLUMNS,
]

HYBRID_REFERENCE_AXIS_FINAL_COLUMNS = [
    *WORKFLOW_AXIS_FINAL_CORE_COLUMNS,
    *WORKFLOW_AXIS_FINAL_LNCRNA_COLUMNS,
]

PAIRED_COHORT_AXIS_FINAL_REQUIRED_COLUMNS = [
    *WORKFLOW_AXIS_FINAL_CORE_COLUMNS,
]

HYBRID_REFERENCE_AXIS_FINAL_REQUIRED_COLUMNS = [
    *WORKFLOW_AXIS_FINAL_CORE_COLUMNS,
]

WORKFLOW_AXIS_FINAL_NUMERIC_COLUMNS = {
    "mRNA_log2FC",
    "miRNA_log2FC",
    "lncRNA_log2FC",
    "circRNA_log2FC",
}


class WorkflowAxisFinalInputError(ValueError):
    pass


class WorkflowAxisFinalPathError(ValueError):
    pass


def validate_safe_name(value: str, field_name: str) -> None:
    value = str(value or "").strip()

    if not value:
        raise WorkflowAxisFinalPathError(
            f"Missing required parameter: {field_name}."
        )

    if "/" in value or "\\" in value or ".." in value:
        raise WorkflowAxisFinalPathError(
            f"Invalid {field_name} parameter."
        )


def get_workflow_task_output_dir(task) -> Path:
    return Path(task.get_output_dir_absolute_path()).resolve()


def get_workflow_axis_final_file_path(task) -> Path:
    task_name = str(task.task_name).strip()

    validate_safe_name(task_name, "task_name")

    output_dir = get_workflow_task_output_dir(task)

    file_path = (
        output_dir / f"{task_name}{WORKFLOW_AXIS_FINAL_FILENAME_SUFFIX}"
    ).resolve()

    if not str(file_path).startswith(str(output_dir)):
        raise WorkflowAxisFinalPathError(
            "Invalid workflow ceRNA axis final file path."
        )

    return file_path


def validate_workflow_axis_final_file(task) -> Path:
    file_path = get_workflow_axis_final_file_path(task)

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Workflow ceRNA axis final file not found: {file_path.name}"
        )

    return file_path


def read_workflow_axis_final_file(
    task,
    required_columns: list[str] | set[str],
) -> tuple[Path, pd.DataFrame]:
    file_path = validate_workflow_axis_final_file(task)

    return read_axis_final_file_by_path(
        file_path=file_path,
        required_columns=required_columns,
    )


def normalize_workflow_axis_final_dataframe(
    df: pd.DataFrame,
    columns: list[str],
    numeric_columns: set[str] | None = None,
) -> pd.DataFrame:
    numeric_columns = numeric_columns or WORKFLOW_AXIS_FINAL_NUMERIC_COLUMNS

    normalized_df = df.copy()

    for col in columns:
        if col not in normalized_df.columns:
            normalized_df[col] = pd.NA

        if col in numeric_columns:
            normalized_df[col] = pd.to_numeric(
                normalized_df[col],
                errors="coerce",
            )
        else:
            normalized_df[col] = normalized_df[col].astype("object")

    normalized_df = normalized_df.replace(
        [float("inf"), float("-inf")],
        pd.NA,
    )

    return normalized_df[columns]


def safe_float_or_none(value):
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
    if value is None or pd.isna(value):
        return ""

    return str(value)


def serialize_workflow_axis_final_dataframe(
    df: pd.DataFrame,
    columns: list[str],
    numeric_columns: set[str] | None = None,
) -> list[dict]:
    numeric_columns = numeric_columns or WORKFLOW_AXIS_FINAL_NUMERIC_COLUMNS

    clean_df = df.replace(
        [float("inf"), float("-inf")],
        pd.NA,
    )

    records = []

    for _, row in clean_df.iterrows():
        item = {}

        for col in columns:
            value = row.get(col)

            if col in numeric_columns:
                item[col] = safe_float_or_none(value)
            else:
                item[col] = safe_str_or_empty(value)

        records.append(item)

    return records


def validate_axis_final_file_path(file_path: Path) -> Path:
    file_path = Path(file_path).resolve()

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"ceRNA axis final file not found: {file_path.name}"
        )

    return file_path


def read_axis_final_file_by_path(
    file_path: Path,
    required_columns: list[str] | set[str],
) -> tuple[Path, pd.DataFrame]:
    file_path = validate_axis_final_file_path(file_path)

    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        raise WorkflowAxisFinalInputError(
            f"Failed to read ceRNA axis final file: {str(e)}"
        )

    required_columns = set(required_columns)
    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise WorkflowAxisFinalInputError(
            "ceRNA axis final file is missing required column(s): "
            f"{', '.join(sorted(missing_columns))}."
        )

    return file_path, df


def build_axis_final_response_from_dataframe(
    *,
    df: pd.DataFrame,
    axis_file_name: str,
    columns: list[str],
    required_columns: list[str] | set[str] | None = None,
    numeric_columns: set[str] | None = None,
    base_response: dict | None = None,
) -> dict:
    numeric_columns = numeric_columns or WORKFLOW_AXIS_FINAL_NUMERIC_COLUMNS
    required_columns = required_columns or WORKFLOW_AXIS_FINAL_CORE_COLUMNS

    visible_columns = resolve_existing_axis_final_columns(
        df=df,
        candidate_columns=columns,
        required_columns=required_columns,
    )

    df = normalize_workflow_axis_final_dataframe(
        df=df,
        columns=visible_columns,
        numeric_columns=numeric_columns,
    )

    raw_count = int(df.shape[0])

    results = serialize_workflow_axis_final_dataframe(
        df=df,
        columns=visible_columns,
        numeric_columns=numeric_columns,
    )

    response_data = {
        "axis_final_file": axis_file_name,
        "count": raw_count,
        "columns": visible_columns,
        "results": results,
    }

    if base_response:
        response_data = {
            **base_response,
            **response_data,
        }

    return response_data


def resolve_existing_axis_final_columns(
    *,
    df: pd.DataFrame,
    candidate_columns: list[str],
    required_columns: list[str] | set[str],
) -> list[str]:
    """
    Return ordered response columns.

    Rules:
    - required columns must exist.
    - optional columns are included only when they exist in the file.
    - column order follows candidate_columns.
    """
    required_columns = set(required_columns)
    existing_columns = set(df.columns)

    missing_required_columns = required_columns - existing_columns

    if missing_required_columns:
        raise WorkflowAxisFinalInputError(
            "ceRNA axis final file is missing required column(s): "
            f"{', '.join(sorted(missing_required_columns))}."
        )

    return [
        col
        for col in candidate_columns
        if col in existing_columns
    ]
