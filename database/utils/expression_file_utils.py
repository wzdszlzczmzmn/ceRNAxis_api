from pathlib import Path

from django.conf import settings


EXPRESSION_TYPES_BY_RNA_TYPE = {
    "mRNA": {"log2count", "log2fpkm", "log2fpkmuq", "log2tpm"},
    "lncRNA": {"log2count", "log2fpkm", "log2fpkmuq", "log2tpm"},
    "miRNA": {"log2rpm"},
    "circRNA": {"count"},
}

EXPRESSION_FILE_FORMATS = {"csv", "parquet"}

DEFAULT_EXPRESSION_FILE_FORMAT = "parquet"

MAX_SELECTED_GENES = 30


class ExpressionPathError(ValueError):
    pass


def validate_dataset(dataset: str) -> None:
    if not dataset:
        raise ExpressionPathError("Missing required parameter: dataset")

    if "/" in dataset or "\\" in dataset or ".." in dataset:
        raise ExpressionPathError("Invalid dataset parameter")


def validate_rna_type(rna_type: str) -> None:
    if rna_type not in EXPRESSION_TYPES_BY_RNA_TYPE:
        raise ExpressionPathError(
            f"Invalid rna_type. Allowed values: {sorted(EXPRESSION_TYPES_BY_RNA_TYPE)}"
        )


def validate_expression_type(rna_type: str, expression_type: str) -> None:
    validate_rna_type(rna_type)

    valid_expression_types = EXPRESSION_TYPES_BY_RNA_TYPE[rna_type]

    if expression_type not in valid_expression_types:
        raise ExpressionPathError(
            f"Invalid expression_type '{expression_type}' for rna_type '{rna_type}'. "
            f"Allowed values: {sorted(valid_expression_types)}"
        )


def validate_expression_file_format(file_format: str) -> None:
    if file_format not in EXPRESSION_FILE_FORMATS:
        raise ExpressionPathError(
            f"Invalid file_format. Allowed values: {sorted(EXPRESSION_FILE_FORMATS)}"
        )


def get_dataset_base_dir() -> Path:
    return Path(settings.DATASET_BASE_DIR).resolve()


def get_expression_filename(
    dataset: str,
    expression_type: str,
    file_format: str,
) -> str:
    return f"{dataset}_{expression_type}_exp.{file_format}"


def get_expression_file_path(
    dataset: str,
    rna_type: str,
    expression_type: str,
    file_format: str = DEFAULT_EXPRESSION_FILE_FORMAT,
) -> Path:
    validate_dataset(dataset)
    validate_expression_type(rna_type, expression_type)
    validate_expression_file_format(file_format)

    base_dir = get_dataset_base_dir()

    file_path = (
        base_dir
        / get_expression_filename(
            dataset=dataset,
            expression_type=expression_type,
            file_format=file_format,
        )
    ).resolve()

    if not str(file_path).startswith(str(base_dir)):
        raise ExpressionPathError("Invalid expression file path")

    return file_path


def get_expression_csv_file_path(
    dataset: str,
    rna_type: str,
    expression_type: str,
) -> Path:
    return get_expression_file_path(
        dataset=dataset,
        rna_type=rna_type,
        expression_type=expression_type,
        file_format="csv",
    )


def get_expression_parquet_file_path(
    dataset: str,
    rna_type: str,
    expression_type: str,
) -> Path:
    return get_expression_file_path(
        dataset=dataset,
        rna_type=rna_type,
        expression_type=expression_type,
        file_format="parquet",
    )


def validate_expression_file(
    dataset: str,
    rna_type: str,
    expression_type: str,
    file_format: str = DEFAULT_EXPRESSION_FILE_FORMAT,
) -> Path:
    file_path = get_expression_file_path(
        dataset=dataset,
        rna_type=rna_type,
        expression_type=expression_type,
        file_format=file_format,
    )

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Expression file not found for dataset '{dataset}', "
            f"rna_type '{rna_type}', expression_type '{expression_type}', "
            f"file_format '{file_format}'."
        )

    return file_path


def validate_expression_csv_file(
    dataset: str,
    rna_type: str,
    expression_type: str,
) -> Path:
    return validate_expression_file(
        dataset=dataset,
        rna_type=rna_type,
        expression_type=expression_type,
        file_format="csv",
    )


def validate_expression_parquet_file(
    dataset: str,
    rna_type: str,
    expression_type: str,
) -> Path:
    return validate_expression_file(
        dataset=dataset,
        rna_type=rna_type,
        expression_type=expression_type,
        file_format="parquet",
    )


def get_available_expression_types(
    dataset: str,
    rna_type: str,
    file_format: str = DEFAULT_EXPRESSION_FILE_FORMAT,
) -> list[str]:
    validate_dataset(dataset)
    validate_rna_type(rna_type)
    validate_expression_file_format(file_format)

    available_types = []

    for expression_type in sorted(EXPRESSION_TYPES_BY_RNA_TYPE[rna_type]):
        file_path = get_expression_file_path(
            dataset=dataset,
            rna_type=rna_type,
            expression_type=expression_type,
            file_format=file_format,
        )

        if file_path.exists() and file_path.is_file():
            available_types.append(expression_type)

    return available_types


def get_available_expression_file_formats(
    dataset: str,
    rna_type: str,
    expression_type: str,
) -> list[str]:
    validate_dataset(dataset)
    validate_expression_type(rna_type, expression_type)

    available_formats = []

    for file_format in sorted(EXPRESSION_FILE_FORMATS):
        file_path = get_expression_file_path(
            dataset=dataset,
            rna_type=rna_type,
            expression_type=expression_type,
            file_format=file_format,
        )

        if file_path.exists() and file_path.is_file():
            available_formats.append(file_format)

    return available_formats
