from pathlib import Path
import csv

PAIRED_COHORT_ALLOWED_FILE_FIELDS = [
    "mrna_file",
    "mirna_file",
    "lncrna_file",
    "meta_file",
]

PAIRED_COHORT_INPUT_FILENAME_MAP = {
    "mrna_file": "mrna.csv",
    "mirna_file": "mirna.csv",
    "lncrna_file": "lncrna.csv",
    "meta_file": "meta.csv",
}

EXPRESSION_REQUIRED_COLUMNS = ["sample_id"]
META_REQUIRED_COLUMNS = ["sample_id", "c_group"]
META_REQUIRED_GROUPS = {"case", "control"}


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
        "output_dir": output_dir
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
