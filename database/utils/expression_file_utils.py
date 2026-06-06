from pathlib import Path

from django.conf import settings


EXPRESSION_FILE_TYPES = {
    "log2count": "{dataset}_log2count_exp.csv",
    "log2fpkm": "{dataset}_log2fpkm_exp.csv",
    "log2fpkmuq": "{dataset}_log2fpkmuq_exp.csv",
    "log2tpm": "{dataset}_log2tpm_exp.csv",
}

MAX_SELECTED_GENES = 30


def get_available_expression_types(dataset: str) -> list[str]:
    available_types = []

    for expression_type, filename_template in EXPRESSION_FILE_TYPES.items():
        filename = filename_template.format(dataset=dataset)
        file_path = settings.DATASET_BASE_DIR / filename

        if file_path.exists() and file_path.is_file():
            available_types.append(expression_type)

    return available_types


def get_expression_file_path(dataset: str, expression_type: str) -> Path:
    if expression_type not in EXPRESSION_FILE_TYPES:
        raise ValueError("Invalid expression_type.")

    filename = EXPRESSION_FILE_TYPES[expression_type].format(dataset=dataset)
    return settings.DATASET_BASE_DIR / filename


def validate_expression_file(dataset: str, expression_type: str) -> Path:
    file_path = get_expression_file_path(dataset, expression_type)

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Expression file not found for dataset '{dataset}' and type '{expression_type}'."
        )

    return file_path
