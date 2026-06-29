import csv
from pathlib import Path
import re

ALLOWED_RNA_TYPES = ["miRNA", "mRNA", "lncRNA", "circRNA"]
MAX_TOTAL_RNA_COUNT = 120

CANCER_TYPE_PATTERN = re.compile(r"^[A-Za-z0-9_]+$")

IMMUNE_FILE_PREFIX = "ImmiRImmiR_"
IMMUNE_FILE_SUFFIX = ".csv"

MIRNA_AXIS_FILE_SUFFIX = "_cerna_axis_miRNA.csv"


class CustomListQueryTaskInputError(ValueError):
    pass


class CustomListQueryPathError(ValueError):
    pass


def normalize_rna_list(value) -> list[str]:
    if not isinstance(value, list):
        raise CustomListQueryTaskInputError("RNA value must be a list.")

    normalized = []
    seen = set()

    for item in value:
        name = str(item).strip()

        if not name:
            continue

        if name in seen:
            continue

        seen.add(name)
        normalized.append(name)

    return normalized


def normalize_rnas(rnas) -> dict[str, list[str]]:
    if not isinstance(rnas, dict):
        raise CustomListQueryTaskInputError("Field 'rnas' must be an object.")

    unknown_keys = set(rnas.keys()) - set(ALLOWED_RNA_TYPES)

    if unknown_keys:
        raise CustomListQueryTaskInputError(
            f"Invalid RNA types: {sorted(unknown_keys)}. "
            f"Allowed values are: {ALLOWED_RNA_TYPES}."
        )

    normalized = {}

    for rna_type in ALLOWED_RNA_TYPES:
        try:
            normalized[rna_type] = normalize_rna_list(
                rnas.get(rna_type, [])
            )
        except CustomListQueryTaskInputError:
            raise CustomListQueryTaskInputError(
                f"Field 'rnas.{rna_type}' must be a list."
            )

    total_count = sum(len(values) for values in normalized.values())

    if total_count == 0:
        raise CustomListQueryTaskInputError(
            "At least one RNA must be provided."
        )

    if total_count > MAX_TOTAL_RNA_COUNT:
        raise CustomListQueryTaskInputError(
            f"At most {MAX_TOTAL_RNA_COUNT} RNAs can be submitted. "
            f"Current count: {total_count}."
        )

    return normalized


def validate_safe_name(value: str, field_name: str) -> None:
    if not value:
        raise CustomListQueryPathError(f"Missing required parameter: {field_name}.")

    if "/" in value or "\\" in value or ".." in value:
        raise CustomListQueryPathError(f"Invalid {field_name} parameter.")


def validate_task_name_for_filename(task_name: str) -> None:
    validate_safe_name(str(task_name).strip(), "task_name")


def get_mirna_axis_filename(task_name: str) -> str:
    task_name = str(task_name).strip()
    validate_task_name_for_filename(task_name)
    return f"{task_name}{MIRNA_AXIS_FILE_SUFFIX}"


def get_task_mirna_axis_file_path(task) -> Path:
    input_dir = get_task_input_dir(task)

    file_path = (
        input_dir / get_mirna_axis_filename(task.task_name)
    ).resolve()

    if not str(file_path).startswith(str(input_dir)):
        raise CustomListQueryPathError("Invalid miRNA axis file path.")

    return file_path


def validate_task_mirna_axis_file(task) -> Path:
    file_path = get_task_mirna_axis_file_path(task)

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"ceRNA axis miRNA file not found: {file_path}"
        )

    return file_path


def validate_cancer_type(cancer_type: str) -> str:
    value = str(cancer_type or "").strip()

    if not value:
        raise CustomListQueryTaskInputError("Missing field: cancer_type.")

    if len(value) > 50:
        raise CustomListQueryTaskInputError(
            "cancer_type must be no more than 50 characters."
        )

    if not CANCER_TYPE_PATTERN.fullmatch(value):
        raise CustomListQueryTaskInputError(
            "Invalid cancer_type. Only letters, numbers and underscores are allowed."
        )

    return value


def get_task_input_dir(task) -> Path:
    return Path(task.get_input_dir_absolute_path()).resolve()


def get_task_output_dir(task) -> Path:
    return Path(task.get_output_dir_absolute_path()).resolve()


def prepare_custom_list_query_workspace(task) -> dict:
    input_dir = get_task_input_dir(task)
    output_dir = get_task_output_dir(task)

    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    return {
        "input_dir": input_dir,
        "output_dir": output_dir,
    }


def get_immune_result_file_name(task) -> str:
    task_name = str(task.task_name).strip()
    validate_task_name_for_filename(task_name)

    return f"{task_name}_map_immune_gene.csv"


def get_immune_result_file_path(task) -> Path:
    filename = get_immune_result_file_name(task)
    output_dir = get_task_output_dir(task)

    return output_dir / filename
