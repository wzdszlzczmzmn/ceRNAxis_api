from pathlib import Path
import csv

import pandas as pd


PAIRED_COHORT_ALLOWED_FILE_FIELDS = [
    "mrna_file",
    "mirna_file",
    "lncrna_file",
    "circrna_file",
    "meta_file",
]

PAIRED_COHORT_INPUT_FILENAME_MAP = {
    "mrna_file": "mrna.csv",
    "mirna_file": "mirna.csv",
    "lncrna_file": "lncrna.csv",
    "circrna_file": "circrna.csv",
    "meta_file": "meta.csv",
}

EXPRESSION_REQUIRED_COLUMNS = ["sample_id"]
META_REQUIRED_COLUMNS = ["sample_id", "c_group"]
META_REQUIRED_GROUPS = {"case", "control"}


# Fixed paired cohort pipeline columns / labels.
PAIRED_COHORT_EXPR_SAMPLE_COL = "sample_id"
PAIRED_COHORT_META_SAMPLE_COL = "sample_id"
PAIRED_COHORT_GROUP_COL = "c_group"
PAIRED_COHORT_CASE_LABEL = "case"


# Output file suffixes.
PAIRED_COHORT_CORRELATION_FILENAME_SUFFIX = "_ceRNA_corr.csv"


# Correlation result file.
PAIRED_COHORT_VALID_CORRELATION_TYPES = [
    "miRNA-mRNA",
    "miRNA-lncRNA",
    "lncRNA-mRNA",
]

PAIRED_COHORT_CORRELATION_REQUIRED_COLUMNS = {
    "gene1",
    "gene2",
    "type",
    "correlation type",
    "correlation",
    "p value",
}

PAIRED_COHORT_TYPE_FILE_FIELD_MAP = {
    "miRNA-mRNA": {
        "gene1_file": "mirna_file",
        "gene2_file": "mrna_file",
    },
    "miRNA-lncRNA": {
        "gene1_file": "mirna_file",
        "gene2_file": "lncrna_file",
    },
    "lncRNA-mRNA": {
        "gene1_file": "lncrna_file",
        "gene2_file": "mrna_file",
    },
}


class PairedCohortTaskInputError(ValueError):
    pass


class PairedCohortTaskPathError(ValueError):
    pass


def validate_safe_name(value: str, field_name: str) -> None:
    if not value:
        raise PairedCohortTaskPathError(f"Missing required parameter: {field_name}.")

    if "/" in value or "\\" in value or ".." in value:
        raise PairedCohortTaskPathError(f"Invalid {field_name} parameter.")


def validate_task_name_for_filename(task_name: str) -> None:
    validate_safe_name(str(task_name).strip(), "task_name")


def get_paired_cohort_task_input_dir(task) -> Path:
    return Path(task.get_input_dir_absolute_path()).resolve()


def get_paired_cohort_task_output_dir(task) -> Path:
    return Path(task.get_output_dir_absolute_path()).resolve()


def get_paired_cohort_input_file_path(task, field_name: str) -> Path:
    if field_name not in PAIRED_COHORT_INPUT_FILENAME_MAP:
        raise PairedCohortTaskPathError(f"Invalid input file field: {field_name}.")

    input_dir = get_paired_cohort_task_input_dir(task)

    file_path = (
        input_dir / PAIRED_COHORT_INPUT_FILENAME_MAP[field_name]
    ).resolve()

    if not str(file_path).startswith(str(input_dir)):
        raise PairedCohortTaskPathError("Invalid paired cohort input file path.")

    return file_path


def prepare_paired_cohort_workspace(task) -> dict:
    input_dir = get_paired_cohort_task_input_dir(task)
    output_dir = get_paired_cohort_task_output_dir(task)

    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    return {
        "input_dir": input_dir,
        "output_dir": output_dir,
    }


def save_paired_cohort_uploaded_input_files(task, files) -> dict:
    input_dir = get_paired_cohort_task_input_dir(task)
    input_dir.mkdir(parents=True, exist_ok=True)

    saved_files = {}

    for field_name in PAIRED_COHORT_ALLOWED_FILE_FIELDS:
        if field_name not in files:
            raise PairedCohortTaskInputError(f"Missing uploaded file: {field_name}.")

        uploaded_file = files[field_name]
        file_path = get_paired_cohort_input_file_path(task, field_name)

        with file_path.open("wb+") as destination:
            for chunk in uploaded_file.chunks():
                destination.write(chunk)

        saved_files[field_name] = file_path.name

    return saved_files


def validate_paired_cohort_input_files(task) -> dict:
    validated_files = {}

    for field_name in PAIRED_COHORT_ALLOWED_FILE_FIELDS:
        file_path = get_paired_cohort_input_file_path(task, field_name)

        if not file_path.exists() or not file_path.is_file():
            raise FileNotFoundError(
                f"Paired cohort input file not found: {file_path}"
            )

        validated_files[field_name] = file_path

    return validated_files


def read_csv_header(file_path: Path) -> list[str]:
    try:
        with file_path.open("r", newline="", encoding="utf-8-sig") as f:
            reader = csv.reader(f)
            header = next(reader, None)
    except UnicodeDecodeError:
        raise PairedCohortTaskInputError(
            f"File must be UTF-8 encoded: {file_path.name}."
        )
    except csv.Error as e:
        raise PairedCohortTaskInputError(
            f"Invalid CSV file: {file_path.name}. {str(e)}"
        )

    if not header:
        raise PairedCohortTaskInputError(
            f"CSV file is empty or missing header: {file_path.name}."
        )

    return [str(col).strip() for col in header]


def validate_required_columns(
    file_path: Path,
    required_columns: list[str],
    file_label: str,
) -> None:
    header = read_csv_header(file_path)
    missing_columns = [
        col for col in required_columns
        if col not in header
    ]

    if missing_columns:
        raise PairedCohortTaskInputError(
            f"{file_label} is missing required column(s): "
            f"{', '.join(missing_columns)}."
        )


def validate_expression_file_columns(file_path: Path, file_label: str) -> None:
    validate_required_columns(
        file_path=file_path,
        required_columns=EXPRESSION_REQUIRED_COLUMNS,
        file_label=file_label,
    )


def validate_meta_file_columns_and_groups(file_path: Path) -> None:
    try:
        with file_path.open("r", newline="", encoding="utf-8-sig") as f:
            reader = csv.reader(f)
            header = next(reader, None)

            if not header:
                raise PairedCohortTaskInputError(
                    "Meta file is empty or missing header."
                )

            normalized_header = [str(col).strip() for col in header]

            missing_columns = [
                col for col in META_REQUIRED_COLUMNS
                if col not in normalized_header
            ]

            if missing_columns:
                raise PairedCohortTaskInputError(
                    "Meta file is missing required column(s): "
                    f"{', '.join(missing_columns)}."
                )

            dict_reader = csv.DictReader(
                f,
                fieldnames=normalized_header,
            )

            observed_required_groups = set()
            row_count = 0

            for row in dict_reader:
                row_count += 1

                group = str(row.get("c_group", "")).strip()

                if group in META_REQUIRED_GROUPS:
                    observed_required_groups.add(group)

    except UnicodeDecodeError:
        raise PairedCohortTaskInputError(
            f"File must be UTF-8 encoded: {file_path.name}."
        )
    except csv.Error as e:
        raise PairedCohortTaskInputError(
            f"Invalid CSV file: {file_path.name}. {str(e)}"
        )

    if row_count == 0:
        raise PairedCohortTaskInputError(
            "Meta file has no data rows."
        )

    missing_groups = META_REQUIRED_GROUPS - observed_required_groups

    if missing_groups:
        raise PairedCohortTaskInputError(
            "Meta file column 'c_group' must contain at least one case "
            "and at least one control sample. "
            f"Missing group(s): {', '.join(sorted(missing_groups))}."
        )


def validate_paired_cohort_file_contents(task) -> dict:
    input_files = validate_paired_cohort_input_files(task)

    validate_expression_file_columns(
        file_path=input_files["mrna_file"],
        file_label="mRNA expression file",
    )

    validate_expression_file_columns(
        file_path=input_files["mirna_file"],
        file_label="miRNA expression file",
    )

    validate_expression_file_columns(
        file_path=input_files["lncrna_file"],
        file_label="lncRNA expression file",
    )

    validate_meta_file_columns_and_groups(
        file_path=input_files["meta_file"],
    )

    return input_files


def get_paired_cohort_correlation_file_path(task) -> Path:
    task_name = str(task.task_name).strip()

    validate_safe_name(task_name, "task_name")

    output_dir = get_paired_cohort_task_output_dir(task)

    file_path = (
        output_dir / f"{task_name}{PAIRED_COHORT_CORRELATION_FILENAME_SUFFIX}"
    ).resolve()

    if not str(file_path).startswith(str(output_dir)):
        raise PairedCohortTaskPathError(
            "Invalid paired cohort correlation file path."
        )

    return file_path


def validate_paired_cohort_correlation_file(task) -> Path:
    file_path = get_paired_cohort_correlation_file_path(task)

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Paired cohort correlation file not found: {file_path.name}"
        )

    return file_path


def read_paired_cohort_correlation_file(task) -> tuple[Path, pd.DataFrame]:
    file_path = validate_paired_cohort_correlation_file(task)

    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        raise PairedCohortTaskInputError(
            f"Failed to read paired cohort correlation file: {str(e)}"
        )

    missing_columns = PAIRED_COHORT_CORRELATION_REQUIRED_COLUMNS - set(df.columns)

    if missing_columns:
        raise PairedCohortTaskInputError(
            "Paired cohort correlation file is missing required column(s): "
            f"{', '.join(sorted(missing_columns))}."
        )

    return file_path, df


def get_paired_cohort_exp_file_fields_by_type(type_value: str) -> dict:
    type_value = str(type_value).strip()

    if type_value not in PAIRED_COHORT_TYPE_FILE_FIELD_MAP:
        raise PairedCohortTaskInputError(
            "Invalid type. Allowed values are: "
            f"{', '.join(PAIRED_COHORT_VALID_CORRELATION_TYPES)}."
        )

    return PAIRED_COHORT_TYPE_FILE_FIELD_MAP[type_value]


def normalize_paired_cohort_correlation_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    normalized_df = df.copy()

    for col in ["gene1", "gene2", "type", "correlation type"]:
        if col in normalized_df.columns:
            normalized_df[col] = normalized_df[col].astype(str).str.strip()

    return normalized_df


def get_paired_cohort_available_correlation_pairs(
    df: pd.DataFrame,
) -> list[dict]:
    normalized_df = normalize_paired_cohort_correlation_dataframe(df)

    normalized_df = normalized_df[
        normalized_df["type"].isin(PAIRED_COHORT_VALID_CORRELATION_TYPES)
    ]

    pair_df = (
        normalized_df[["gene1", "gene2", "type"]]
        .dropna()
        .drop_duplicates()
        .sort_values(["type", "gene1", "gene2"])
    )

    return [
        {
            "gene1": row["gene1"],
            "gene2": row["gene2"],
            "type": row["type"],
        }
        for _, row in pair_df.iterrows()
    ]


def extract_paired_cohort_correlation_stats(
    cor_df: pd.DataFrame,
) -> dict:
    stats = {
        "pearson_r": None,
        "pearson_p": None,
        "spearman_r": None,
        "spearman_p": None,
        "kendall_tau": None,
        "kendall_p": None,
    }

    normalized_df = normalize_paired_cohort_correlation_dataframe(cor_df)

    pearson_df = normalized_df[
        normalized_df["correlation type"] == "Pearson Correlation"
    ]

    if not pearson_df.empty:
        row = pearson_df.iloc[0]
        stats["pearson_r"] = safe_float_or_none(row.get("correlation"))
        stats["pearson_p"] = safe_float_or_none(row.get("p value"))

    spearman_df = normalized_df[
        normalized_df["correlation type"] == "Spearman Correlation"
    ]

    if not spearman_df.empty:
        row = spearman_df.iloc[0]
        stats["spearman_r"] = safe_float_or_none(row.get("correlation"))
        stats["spearman_p"] = safe_float_or_none(row.get("p value"))

    kendall_df = normalized_df[
        normalized_df["correlation type"] == "Kendall's tau"
    ]

    if not kendall_df.empty:
        row = kendall_df.iloc[0]
        stats["kendall_tau"] = safe_float_or_none(row.get("correlation"))
        stats["kendall_p"] = safe_float_or_none(row.get("p value"))

    return stats


def safe_float_or_none(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None

    if pd.isna(result):
        return None

    return result
