from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq

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

TIMEDB_RNA_TYPE = "mRNA"
TIMEDB_EXPRESSION_TYPE = "exp"

TISCH2_RNA_TYPE = "mRNA"
TISCH2_EXPRESSION_TYPE = "exp"
TISCH2_EXPRESSION_FILE_FORMAT = "parquet"

SCTML_RNA_TYPE = "mRNA"
SCTML_EXPRESSION_TYPE = "exp"
SCTML_EXPRESSION_FILE_FORMAT = "parquet"
SCTML_EXPRESSION_MODE = "scTML"

ALIQUOT_EXPRESSION_VALUE_TYPES_BY_RNA_TYPE = {
    "mRNA": {"count", "tpm", "fpkm", "fpkmuq"},
    "lncRNA": {"count", "tpm", "fpkm", "fpkmuq"},
    "miRNA": {"count", "rpm"},
    "isoform": {"count", "rpm"},
}

ALIQUOT_EXPRESSION_FILE_FORMAT = "csv"
ALIQUOT_SAMPLE_LEVEL = "aliquot"


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
    return Path(settings.TCGA_DATASET_BASE_DIR).resolve()


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


def get_timedb_dataset_base_dir() -> Path:
    return Path(settings.TIMEDB_DATASET_BASE_DIR).resolve()


def validate_timedb_rna_type(rna_type: str) -> None:
    if rna_type != TIMEDB_RNA_TYPE:
        raise ExpressionPathError(
            f"Invalid TIMEDB rna_type. Allowed value: '{TIMEDB_RNA_TYPE}'"
        )


def get_timedb_expression_filename(
    dataset: str,
    file_format: str,
) -> str:
    return f"{dataset}_exp.{file_format}"


def get_timedb_expression_file_path(
    dataset: str,
    rna_type: str,
    file_format: str = DEFAULT_EXPRESSION_FILE_FORMAT,
) -> Path:
    validate_dataset(dataset)
    validate_timedb_rna_type(rna_type)
    validate_expression_file_format(file_format)

    base_dir = get_timedb_dataset_base_dir()

    file_path = (
        base_dir
        / get_timedb_expression_filename(
            dataset=dataset,
            file_format=file_format,
        )
    ).resolve()

    if not str(file_path).startswith(str(base_dir)):
        raise ExpressionPathError("Invalid TIMEDB expression file path")

    return file_path


def validate_timedb_expression_file(
    dataset: str,
    rna_type: str,
    file_format: str = DEFAULT_EXPRESSION_FILE_FORMAT,
) -> Path:
    file_path = get_timedb_expression_file_path(
        dataset=dataset,
        rna_type=rna_type,
        file_format=file_format,
    )

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"TIMEDB expression file not found for dataset '{dataset}', "
            f"rna_type '{rna_type}', file_format '{file_format}'."
        )

    return file_path


def get_available_timedb_expression_file_formats(
    dataset: str,
    rna_type: str,
) -> list[str]:
    validate_dataset(dataset)
    validate_timedb_rna_type(rna_type)

    available_formats = []

    for file_format in sorted(EXPRESSION_FILE_FORMATS):
        file_path = get_timedb_expression_file_path(
            dataset=dataset,
            rna_type=rna_type,
            file_format=file_format,
        )

        if file_path.exists() and file_path.is_file():
            available_formats.append(file_format)

    return available_formats


def get_default_timedb_expression_file_format(
    dataset: str,
    rna_type: str,
) -> str:
    available_formats = get_available_timedb_expression_file_formats(
        dataset=dataset,
        rna_type=rna_type,
    )

    if DEFAULT_EXPRESSION_FILE_FORMAT in available_formats:
        return DEFAULT_EXPRESSION_FILE_FORMAT

    if available_formats:
        return available_formats[0]

    return DEFAULT_EXPRESSION_FILE_FORMAT


def get_available_timedb_expression_types(
    dataset: str,
    rna_type: str,
    file_format: str | None = None,
) -> list[str]:
    validate_dataset(dataset)
    validate_timedb_rna_type(rna_type)

    if file_format is None:
        file_format = get_default_timedb_expression_file_format(
            dataset=dataset,
            rna_type=rna_type,
        )

    validate_expression_file_format(file_format)

    file_path = get_timedb_expression_file_path(
        dataset=dataset,
        rna_type=rna_type,
        file_format=file_format,
    )

    if file_path.exists() and file_path.is_file():
        return [TIMEDB_EXPRESSION_TYPE]

    return []


def get_expression_mode_from_metadata(metadata) -> str:
    obs_type = str(metadata.obs_type or "").strip().lower()
    programme = str(metadata.programme or "").strip().upper()

    if obs_type == "cell":
        return "tisch2"

    if obs_type == "spot":
        return "scTML"

    if programme == "TCGA":
        return "tcga"

    return "timedb"


def resolve_dataset_expression_file(
    dataset: str,
    rna_type: str,
    expression_type: str,
    expression_mode: str,
):
    if expression_mode == "tcga":
        validate_expression_type(
            rna_type=rna_type,
            expression_type=expression_type,
        )

        file_format = DEFAULT_EXPRESSION_FILE_FORMAT

        file_path = validate_expression_file(
            dataset=dataset,
            rna_type=rna_type,
            expression_type=expression_type,
            file_format=file_format,
        )

        return file_path, file_format

    if expression_mode == "timedb":
        if expression_type != TIMEDB_EXPRESSION_TYPE:
            raise ExpressionPathError(
                f"Invalid TIMEDB expression_type '{expression_type}'. "
                f"Allowed value: '{TIMEDB_EXPRESSION_TYPE}'."
            )

        file_format = get_default_timedb_expression_file_format(
            dataset=dataset,
            rna_type=rna_type,
        )

        file_path = validate_timedb_expression_file(
            dataset=dataset,
            rna_type=rna_type,
            file_format=file_format,
        )

        return file_path, file_format

    if expression_mode == "scTML":
        if expression_type != SCTML_EXPRESSION_TYPE:
            raise ExpressionPathError(
                f"Invalid scTML expression_type '{expression_type}'. "
                f"Allowed value: '{SCTML_EXPRESSION_TYPE}'."
            )

        file_format = SCTML_EXPRESSION_FILE_FORMAT

        file_path = validate_sctml_expression_file(
            dataset=dataset,
            rna_type=rna_type,
        )

        return file_path, file_format

    if expression_mode == "tisch2":
        if expression_type != TISCH2_EXPRESSION_TYPE:
            raise ExpressionPathError(
                f"Invalid TISCH2 expression_type '{expression_type}'. "
                f"Allowed value: '{TISCH2_EXPRESSION_TYPE}'."
            )

        file_format = TISCH2_EXPRESSION_FILE_FORMAT

        file_path = validate_tisch2_expression_file(
            dataset=dataset,
            rna_type=rna_type,
        )

        return file_path, file_format

    raise ExpressionPathError(
        f"Invalid expression_mode '{expression_mode}'."
    )


def read_expression_columns(
    file_path: Path,
    file_format: str,
) -> list[str]:
    if file_format == "parquet":
        schema = pq.read_schema(file_path)
        return schema.names

    if file_format == "csv":
        df_header = pd.read_csv(
            file_path,
            nrows=0,
            low_memory=False,
        )
        return df_header.columns.tolist()

    raise ExpressionPathError(
        f"Unsupported expression file_format '{file_format}'."
    )


def read_expression_data(
    file_path: Path,
    file_format: str,
    usecols: list[str],
) -> pd.DataFrame:
    if file_format == "parquet":
        return pd.read_parquet(
            file_path,
            columns=usecols,
            engine="pyarrow",
        )

    if file_format == "csv":
        return pd.read_csv(
            file_path,
            usecols=usecols,
            low_memory=False,
        )

    raise ExpressionPathError(
        f"Unsupported expression file_format '{file_format}'."
    )


def get_tisch2_dataset_base_dir() -> Path:
    return Path(settings.TISCH2_DATASET_BASE_DIR).resolve()


def validate_tisch2_rna_type(rna_type: str) -> None:
    if rna_type != TISCH2_RNA_TYPE:
        raise ExpressionPathError(
            f"Invalid TISCH2 rna_type. Allowed value: '{TISCH2_RNA_TYPE}'"
        )


def get_tisch2_expression_filename(dataset: str) -> str:
    return f"{dataset}_exp.{TISCH2_EXPRESSION_FILE_FORMAT}"


def get_tisch2_expression_file_path(
    dataset: str,
    rna_type: str,
) -> Path:
    validate_dataset(dataset)
    validate_tisch2_rna_type(rna_type)

    base_dir = get_tisch2_dataset_base_dir()

    file_path = (
        base_dir
        / get_tisch2_expression_filename(dataset=dataset)
    ).resolve()

    if not str(file_path).startswith(str(base_dir)):
        raise ExpressionPathError("Invalid TISCH2 expression file path")

    return file_path


def validate_tisch2_expression_file(
    dataset: str,
    rna_type: str,
) -> Path:
    file_path = get_tisch2_expression_file_path(
        dataset=dataset,
        rna_type=rna_type,
    )

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"TISCH2 expression file not found for dataset '{dataset}', "
            f"rna_type '{rna_type}', file_format '{TISCH2_EXPRESSION_FILE_FORMAT}'."
        )

    return file_path


def get_available_tisch2_expression_types(
    dataset: str,
    rna_type: str,
) -> list[str]:
    validate_dataset(dataset)
    validate_tisch2_rna_type(rna_type)

    file_path = get_tisch2_expression_file_path(
        dataset=dataset,
        rna_type=rna_type,
    )

    if file_path.exists() and file_path.is_file():
        return [TISCH2_EXPRESSION_TYPE]

    return []


def get_sctml_dataset_base_dir() -> Path:
    return Path(settings.SCTML_DATASET_BASE_DIR).resolve()


def validate_sctml_rna_type(rna_type: str) -> None:
    if rna_type != SCTML_RNA_TYPE:
        raise ExpressionPathError(
            f"Invalid scTML rna_type. Allowed value: '{SCTML_RNA_TYPE}'"
        )


def get_sctml_expression_filename(dataset: str) -> str:
    return f"{dataset}_exp.{SCTML_EXPRESSION_FILE_FORMAT}"


def get_sctml_expression_file_path(
    dataset: str,
    rna_type: str,
) -> Path:
    validate_dataset(dataset)
    validate_sctml_rna_type(rna_type)

    base_dir = get_sctml_dataset_base_dir()

    file_path = (
        base_dir
        / get_sctml_expression_filename(dataset=dataset)
    ).resolve()

    if not str(file_path).startswith(str(base_dir)):
        raise ExpressionPathError("Invalid scTML expression file path")

    return file_path


def validate_sctml_expression_file(
    dataset: str,
    rna_type: str,
) -> Path:
    file_path = get_sctml_expression_file_path(
        dataset=dataset,
        rna_type=rna_type,
    )

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"scTML expression file not found for dataset '{dataset}', "
            f"rna_type '{rna_type}', file_format '{SCTML_EXPRESSION_FILE_FORMAT}'."
        )

    return file_path


def get_available_sctml_expression_types(
    dataset: str,
    rna_type: str,
) -> list[str]:
    validate_dataset(dataset)
    validate_sctml_rna_type(rna_type)

    file_path = get_sctml_expression_file_path(
        dataset=dataset,
        rna_type=rna_type,
    )

    if file_path.exists() and file_path.is_file():
        return [SCTML_EXPRESSION_TYPE]

    return []


def validate_aliquot_rna_type(rna_type: str) -> None:
    if rna_type not in ALIQUOT_EXPRESSION_VALUE_TYPES_BY_RNA_TYPE:
        raise ExpressionPathError(
            f"Invalid aliquot rna_type. "
            f"Allowed values: {sorted(ALIQUOT_EXPRESSION_VALUE_TYPES_BY_RNA_TYPE)}"
        )


def validate_aliquot_value_type(rna_type: str, value_type: str) -> None:
    validate_aliquot_rna_type(rna_type)

    valid_value_types = ALIQUOT_EXPRESSION_VALUE_TYPES_BY_RNA_TYPE[rna_type]

    if value_type not in valid_value_types:
        raise ExpressionPathError(
            f"Invalid aliquot value_type '{value_type}' for rna_type '{rna_type}'. "
            f"Allowed values: {sorted(valid_value_types)}"
        )


def get_aliquot_expression_filename(
    dataset: str,
    value_type: str,
) -> str:
    return (
        f"{dataset}_{value_type}_"
        f"{ALIQUOT_SAMPLE_LEVEL}_exp.{ALIQUOT_EXPRESSION_FILE_FORMAT}"
    )


def get_aliquot_isoform_info_filename(dataset: str) -> str:
    return (
        f"{dataset}_info_"
        f"{ALIQUOT_SAMPLE_LEVEL}.{ALIQUOT_EXPRESSION_FILE_FORMAT}"
    )


def get_aliquot_expression_file_path(
    dataset: str,
    rna_type: str,
    value_type: str,
) -> Path:
    validate_dataset(dataset)
    validate_aliquot_value_type(rna_type, value_type)

    base_dir = get_dataset_base_dir()

    file_path = (
        base_dir
        / get_aliquot_expression_filename(
            dataset=dataset,
            value_type=value_type,
        )
    ).resolve()

    if not str(file_path).startswith(str(base_dir)):
        raise ExpressionPathError("Invalid aliquot expression file path")

    return file_path


def get_aliquot_isoform_info_file_path(dataset: str) -> Path:
    validate_dataset(dataset)

    base_dir = get_dataset_base_dir()

    file_path = (
        base_dir
        / get_aliquot_isoform_info_filename(dataset=dataset)
    ).resolve()

    if not str(file_path).startswith(str(base_dir)):
        raise ExpressionPathError("Invalid aliquot isoform info file path")

    return file_path


def get_forced_isoform_dataset_from_mirna_dataset(dataset: str) -> str:
    suffix = "_miRNA"

    if not dataset.endswith(suffix):
        raise ExpressionPathError(
            f"Cannot infer isoform dataset from non-miRNA dataset '{dataset}'."
        )

    return f"{dataset[:-len(suffix)]}_isoform"


def get_available_aliquot_expression_files(
    dataset: str,
    rna_type: str,
) -> list[dict]:
    validate_dataset(dataset)
    validate_aliquot_rna_type(rna_type)

    results = []

    for value_type in sorted(ALIQUOT_EXPRESSION_VALUE_TYPES_BY_RNA_TYPE[rna_type]):
        file_path = get_aliquot_expression_file_path(
            dataset=dataset,
            rna_type=rna_type,
            value_type=value_type,
        )

        if file_path.exists() and file_path.is_file():
            results.append(
                {
                    "file_id": f"{dataset}__{value_type}__aliquot_csv",
                    "dataset": dataset,
                    "filename": file_path.name,
                    "source_rna_type": rna_type,
                    "value_type": value_type,
                    "file_type": "expression",
                    "file_format": ALIQUOT_EXPRESSION_FILE_FORMAT,
                    "sample_level": ALIQUOT_SAMPLE_LEVEL,
                }
            )

    if rna_type == "isoform":
        file_path = get_aliquot_isoform_info_file_path(dataset=dataset)

        if file_path.exists() and file_path.is_file():
            results.append(
                {
                    "file_id": f"{dataset}__info__aliquot_csv",
                    "dataset": dataset,
                    "filename": file_path.name,
                    "source_rna_type": "isoform",
                    "value_type": None,
                    "file_type": "annotation",
                    "file_format": ALIQUOT_EXPRESSION_FILE_FORMAT,
                    "sample_level": ALIQUOT_SAMPLE_LEVEL,
                }
            )

    return results
