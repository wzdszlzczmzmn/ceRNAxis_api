from pathlib import Path
import csv
import math

import pandas as pd

from analysis.utils.workflow_detail_utils.deg_pathway_utils import DEG_PATHWAY_REQUIRED_COLUMNS, \
    validate_deg_pathway_dataframe_columns, DEGPathwayInputError, build_deg_pathway_data_from_dataframe
from analysis.utils.workflow_detail_utils.survival_km_utils import SURVIVAL_REQUIRED_COLUMNS, DEFAULT_SURVIVAL_GROUPS, \
    validate_survival_dataframe_columns, SurvivalKMInputError, build_survival_km_data_from_dataframe

PAIRED_COHORT_REQUIRED_FILE_FIELDS = [
    "mrna_file",
    "mirna_file",
    "meta_file",
]

PAIRED_COHORT_OPTIONAL_CERNA_FILE_FIELDS = [
    "lncrna_file",
    "circrna_file",
]

PAIRED_COHORT_ALLOWED_FILE_FIELDS = [
    *PAIRED_COHORT_REQUIRED_FILE_FIELDS,
    *PAIRED_COHORT_OPTIONAL_CERNA_FILE_FIELDS,
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


PAIRED_COHORT_DEG_RNA_TYPES = [
    "mRNA",
    "miRNA",
    "lncRNA",
    "circRNA",
]


# Background result file.
PAIRED_COHORT_BACKGROUND_FILENAME_SUFFIX = "_ceRNA_background.csv"

PAIRED_COHORT_VALID_BACKGROUND_TYPES = [
    "miRNA-mRNA",
    "miRNA-lncRNA",
    "miRNA-circRNA",
]

PAIRED_COHORT_BACKGROUND_REQUIRED_COLUMNS = {
    "miRNA",
    "ceRNA",
    "species",
    "database",
    "type",
    "miRNA_log2FC",
    "ceRNA_log2FC",
}


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


# ceRNA axis final result file.
PAIRED_COHORT_AXIS_FINAL_FILENAME_SUFFIX = "_ceRNA_axis_final.csv"

PAIRED_COHORT_AXIS_FINAL_COLUMNS = [
    "axis_id",
    "axis_regulation",
    "axis_type",
    "mRNA",
    "mRNA_log2FC",
    "mRNA_regulation",
    "miRNA",
    "miRNA_log2FC",
    "miRNA_regulation",
    "lncRNA",
    "lncRNA_log2FC",
    "lncRNA_regulation",
    "circRNA",
    "circRNA_log2FC",
    "circRNA_regulation",
]

PAIRED_COHORT_AXIS_FINAL_REQUIRED_COLUMNS = set(
    PAIRED_COHORT_AXIS_FINAL_COLUMNS
)


# Survival analysis result file.
PAIRED_COHORT_SURVIVAL_FILENAME_SUFFIX = "_survival_analysis.csv"

PAIRED_COHORT_SURVIVAL_REQUIRED_COLUMNS = SURVIVAL_REQUIRED_COLUMNS

PAIRED_COHORT_SURVIVAL_GROUPS = DEFAULT_SURVIVAL_GROUPS


# DEG pathway / GSEA result file.
PAIRED_COHORT_MRNA_GSEA_FILENAME_SUFFIX = "_mRNA_gsea.csv"

PAIRED_COHORT_MRNA_GSEA_REQUIRED_COLUMNS = DEG_PATHWAY_REQUIRED_COLUMNS


# CMap result file.
PAIRED_COHORT_CMAP_FILENAME_SUFFIX = "_CMap.csv"

PAIRED_COHORT_CMAP_REQUIRED_COLUMNS = set()


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

    missing_required_files = [
        field_name
        for field_name in PAIRED_COHORT_REQUIRED_FILE_FIELDS
        if field_name not in files
    ]

    if missing_required_files:
        raise PairedCohortTaskInputError(
            "Missing uploaded file(s): "
            f"{', '.join(missing_required_files)}."
        )

    has_lncrna_file = "lncrna_file" in files
    has_circrna_file = "circrna_file" in files

    if not has_lncrna_file and not has_circrna_file:
        raise PairedCohortTaskInputError(
            "At least one of lncrna_file or circrna_file is required."
        )

    for field_name in PAIRED_COHORT_ALLOWED_FILE_FIELDS:
        if field_name not in files:
            saved_files[field_name] = ""
            continue

        uploaded_file = files[field_name]
        file_path = get_paired_cohort_input_file_path(task, field_name)

        with file_path.open("wb+") as destination:
            for chunk in uploaded_file.chunks():
                destination.write(chunk)

        saved_files[field_name] = file_path.name

    return saved_files


def validate_paired_cohort_input_files(task) -> dict:
    validated_files = {}

    for field_name in PAIRED_COHORT_REQUIRED_FILE_FIELDS:
        file_path = get_paired_cohort_input_file_path(task, field_name)

        if not file_path.exists() or not file_path.is_file():
            raise FileNotFoundError(
                f"Paired cohort input file not found: {file_path}"
            )

        validated_files[field_name] = file_path

    optional_existing_files = []

    for field_name in PAIRED_COHORT_OPTIONAL_CERNA_FILE_FIELDS:
        file_path = get_paired_cohort_input_file_path(task, field_name)

        if file_path.exists() and file_path.is_file():
            validated_files[field_name] = file_path
            optional_existing_files.append(field_name)
        else:
            validated_files[field_name] = None

    if not optional_existing_files:
        raise FileNotFoundError(
            "At least one paired cohort ceRNA input file is required: "
            "lncrna_file or circrna_file."
        )

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

    if input_files.get("lncrna_file") is not None:
        validate_expression_file_columns(
            file_path=input_files["lncrna_file"],
            file_label="lncRNA expression file",
        )

    if input_files.get("circrna_file") is not None:
        validate_expression_file_columns(
            file_path=input_files["circrna_file"],
            file_label="circRNA expression file",
        )

    validate_meta_file_columns_and_groups(
        file_path=input_files["meta_file"],
    )

    return input_files


def get_paired_cohort_deg_file_path(task, rna_type: str) -> Path:
    rna_type = str(rna_type).strip()

    if rna_type not in PAIRED_COHORT_DEG_RNA_TYPES:
        raise PairedCohortTaskPathError(
            "Invalid DEG rna_type. Allowed values are: "
            f"{', '.join(PAIRED_COHORT_DEG_RNA_TYPES)}."
        )

    task_name = str(task.task_name).strip()
    deg_method = str(task.deg_method).strip()

    validate_safe_name(task_name, "task_name")
    validate_safe_name(deg_method, "deg_method")

    output_dir = get_paired_cohort_task_output_dir(task)

    file_path = (
        output_dir / f"{task_name}_{deg_method}_{rna_type}.csv"
    ).resolve()

    if not str(file_path).startswith(str(output_dir)):
        raise PairedCohortTaskPathError(
            "Invalid paired cohort DEG file path."
        )

    return file_path


def get_available_paired_cohort_deg_rna_types(task) -> list[str]:
    output_dir = get_paired_cohort_task_output_dir(task)

    if not output_dir.exists() or not output_dir.is_dir():
        return []

    available_rna_types = []

    for rna_type in PAIRED_COHORT_DEG_RNA_TYPES:
        file_path = get_paired_cohort_deg_file_path(
            task=task,
            rna_type=rna_type,
        )

        if file_path.exists() and file_path.is_file():
            available_rna_types.append(rna_type)

    return available_rna_types


def get_paired_cohort_background_file_path(task) -> Path:
    task_name = str(task.task_name).strip()

    validate_safe_name(task_name, "task_name")

    output_dir = get_paired_cohort_task_output_dir(task)

    file_path = (
        output_dir / f"{task_name}{PAIRED_COHORT_BACKGROUND_FILENAME_SUFFIX}"
    ).resolve()

    if not str(file_path).startswith(str(output_dir)):
        raise PairedCohortTaskPathError(
            "Invalid paired cohort ceRNA background file path."
        )

    return file_path


def validate_paired_cohort_background_file(task) -> Path:
    file_path = get_paired_cohort_background_file_path(task)

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Paired cohort ceRNA background file not found: {file_path.name}"
        )

    return file_path


def read_paired_cohort_background_file(task) -> tuple[Path, pd.DataFrame]:
    file_path = validate_paired_cohort_background_file(task)

    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        raise PairedCohortTaskInputError(
            f"Failed to read paired cohort ceRNA background file: {str(e)}"
        )

    missing_columns = PAIRED_COHORT_BACKGROUND_REQUIRED_COLUMNS - set(df.columns)

    if missing_columns:
        raise PairedCohortTaskInputError(
            "Paired cohort ceRNA background file is missing required column(s): "
            f"{', '.join(sorted(missing_columns))}."
        )

    return file_path, df


def normalize_paired_cohort_background_dataframe(
    df: pd.DataFrame,
) -> pd.DataFrame:
    normalized_df = df.copy()

    for col in ["miRNA", "ceRNA", "species", "database", "type"]:
        if col in normalized_df.columns:
            normalized_df[col] = normalized_df[col].astype(str).str.strip()

    return normalized_df


def get_paired_cohort_available_background_types(
    df: pd.DataFrame,
) -> list[str]:
    normalized_df = normalize_paired_cohort_background_dataframe(df)

    observed_types = set(
        normalized_df["type"]
        .dropna()
        .astype(str)
        .str.strip()
        .tolist()
    )

    return [
        type_value
        for type_value in PAIRED_COHORT_VALID_BACKGROUND_TYPES
        if type_value in observed_types
    ]


def get_available_paired_cohort_background_types(task) -> list[str]:
    try:
        _, df = read_paired_cohort_background_file(task)
    except (FileNotFoundError, PairedCohortTaskInputError):
        return []

    return get_paired_cohort_available_background_types(df)


def get_paired_cohort_survival_file_path(task) -> Path:
    task_name = str(task.task_name).strip()

    validate_safe_name(task_name, "task_name")

    output_dir = get_paired_cohort_task_output_dir(task)

    file_path = (
        output_dir / f"{task_name}{PAIRED_COHORT_SURVIVAL_FILENAME_SUFFIX}"
    ).resolve()

    if not str(file_path).startswith(str(output_dir)):
        raise PairedCohortTaskPathError(
            "Invalid paired cohort survival analysis file path."
        )

    return file_path


def validate_paired_cohort_survival_file(task) -> Path:
    file_path = get_paired_cohort_survival_file_path(task)

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Paired cohort survival analysis file not found: {file_path.name}"
        )

    return file_path


def read_paired_cohort_survival_file(task) -> tuple[Path, pd.DataFrame]:
    file_path = validate_paired_cohort_survival_file(task)

    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        raise PairedCohortTaskInputError(
            f"Failed to read paired cohort survival analysis file: {str(e)}"
        )

    try:
        validate_survival_dataframe_columns(
            df=df,
            required_columns=PAIRED_COHORT_SURVIVAL_REQUIRED_COLUMNS,
        )
    except SurvivalKMInputError as e:
        raise PairedCohortTaskInputError(
            "Paired cohort survival analysis file is invalid. "
            f"{str(e)}"
        )

    return file_path, df


def build_paired_cohort_survival_km_data(
    task,
    title: str = "ceRNA axis-based survival analysis",
) -> dict:
    survival_file, df = read_paired_cohort_survival_file(task)

    try:
        return build_survival_km_data_from_dataframe(
            task=task,
            survival_file_name=survival_file.name,
            df=df,
            title=title,
            valid_groups=PAIRED_COHORT_SURVIVAL_GROUPS,
        )
    except SurvivalKMInputError as e:
        raise PairedCohortTaskInputError(str(e))


def get_paired_cohort_mrna_gsea_file_path(task) -> Path:
    task_name = str(task.task_name).strip()

    validate_safe_name(task_name, "task_name")

    output_dir = get_paired_cohort_task_output_dir(task)

    file_path = (
        output_dir / f"{task_name}{PAIRED_COHORT_MRNA_GSEA_FILENAME_SUFFIX}"
    ).resolve()

    if not str(file_path).startswith(str(output_dir)):
        raise PairedCohortTaskPathError(
            "Invalid paired cohort mRNA GSEA result file path."
        )

    return file_path


def validate_paired_cohort_mrna_gsea_file(task) -> Path:
    file_path = get_paired_cohort_mrna_gsea_file_path(task)

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Paired cohort mRNA GSEA result file not found: {file_path.name}"
        )

    return file_path


def read_paired_cohort_mrna_gsea_file(task) -> tuple[Path, pd.DataFrame]:
    file_path = validate_paired_cohort_mrna_gsea_file(task)

    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        raise PairedCohortTaskInputError(
            f"Failed to read paired cohort mRNA GSEA result file: {str(e)}"
        )

    try:
        validate_deg_pathway_dataframe_columns(
            df=df,
            required_columns=PAIRED_COHORT_MRNA_GSEA_REQUIRED_COLUMNS,
        )
    except DEGPathwayInputError as e:
        raise PairedCohortTaskInputError(
            "Paired cohort mRNA GSEA result file is invalid. "
            f"{str(e)}"
        )

    return file_path, df


def build_paired_cohort_deg_pathway_data(
    task,
    title: str = "DEG Pathway Enrichment",
) -> dict:
    gsea_file, df = read_paired_cohort_mrna_gsea_file(task)

    try:
        return build_deg_pathway_data_from_dataframe(
            task=task,
            gsea_file_name=gsea_file.name,
            df=df,
            title=title,
        )
    except DEGPathwayInputError as e:
        raise PairedCohortTaskInputError(str(e))
