from pathlib import Path

from django.conf import settings


DEG_EXPRESSION_TYPES_BY_RNA_TYPE = {
    "mRNA": {"log2count", "log2fpkm", "log2fpkmuq", "log2tpm"},
    "lncRNA": {"log2count", "log2fpkm", "log2fpkmuq", "log2tpm"},
    "miRNA": {"log2rpm"},
    "circRNA": {"count"},
}

TIMEDB_RNA_TYPE = "mRNA"
TIMEDB_EXPRESSION_TYPE = "exp"
TCGA_PROGRAMME = "TCGA"

TISCH2_RNA_TYPE = "mRNA"
TISCH2_EXPRESSION_TYPE = "exp"
CELL_OBS_TYPE = "cell"
SPOT_OBS_TYPE = "spot"


def is_tcga_programme(programme: str) -> bool:
    return str(programme).upper() == TCGA_PROGRAMME


def is_tcga_dataset_metadata(metadata) -> bool:
    return is_tcga_programme(metadata.programme)


def is_cell_dataset_metadata(metadata) -> bool:
    return str(metadata.obs_type or "").strip().lower() == CELL_OBS_TYPE


def is_spot_dataset_metadata(metadata) -> bool:
    return str(metadata.obs_type or "").strip().lower() == SPOT_OBS_TYPE


class DEGPathError(ValueError):
    pass


def validate_dataset(dataset: str) -> None:
    if not dataset:
        raise DEGPathError("Missing required parameter: dataset")

    if "/" in dataset or "\\" in dataset or ".." in dataset:
        raise DEGPathError("Invalid dataset parameter")


def validate_rna_type(rna_type: str) -> None:
    if rna_type not in DEG_EXPRESSION_TYPES_BY_RNA_TYPE:
        raise DEGPathError(
            f"Invalid rna_type. Allowed values: {sorted(DEG_EXPRESSION_TYPES_BY_RNA_TYPE)}"
        )


def validate_expression_type(rna_type: str, expression_type: str) -> None:
    validate_rna_type(rna_type)

    valid_expression_types = DEG_EXPRESSION_TYPES_BY_RNA_TYPE[rna_type]

    if expression_type not in valid_expression_types:
        raise DEGPathError(
            f"Invalid expression_type '{expression_type}' for rna_type '{rna_type}'. "
            f"Allowed values: {sorted(valid_expression_types)}"
        )


def get_dataset_base_dir() -> Path:
    return Path(settings.TCGA_DATASET_BASE_DIR).resolve()


def get_deg_filename(dataset: str, expression_type: str) -> str:
    return f"{dataset}_{expression_type}_deg.csv"


def get_deg_file_path(
    dataset: str,
    rna_type: str,
    expression_type: str,
) -> Path:
    validate_dataset(dataset)
    validate_expression_type(
        rna_type=rna_type,
        expression_type=expression_type,
    )

    base_dir = get_dataset_base_dir()

    file_path = (
        base_dir
        / get_deg_filename(
            dataset=dataset,
            expression_type=expression_type,
        )
    ).resolve()

    if not str(file_path).startswith(str(base_dir)):
        raise DEGPathError("Invalid DEG file path")

    return file_path


def validate_deg_file(
    dataset: str,
    rna_type: str,
    expression_type: str,
) -> Path:
    file_path = get_deg_file_path(
        dataset=dataset,
        rna_type=rna_type,
        expression_type=expression_type,
    )

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"DEG file not found for dataset '{dataset}', "
            f"rna_type '{rna_type}', expression_type '{expression_type}'."
        )

    return file_path


def get_available_deg_expression_types(
    dataset: str,
    rna_type: str,
) -> list[str]:
    validate_dataset(dataset)
    validate_rna_type(rna_type)

    available_expression_types = []

    for expression_type in sorted(DEG_EXPRESSION_TYPES_BY_RNA_TYPE[rna_type]):
        file_path = get_deg_file_path(
            dataset=dataset,
            rna_type=rna_type,
            expression_type=expression_type,
        )

        if file_path.exists() and file_path.is_file():
            available_expression_types.append(expression_type)

    return available_expression_types


def get_timedb_dataset_base_dir() -> Path:
    return Path(settings.TIMEDB_DATASET_BASE_DIR).resolve()


def validate_timedb_rna_type(rna_type: str) -> None:
    if rna_type != TIMEDB_RNA_TYPE:
        raise DEGPathError(
            f"Invalid TIMEDB rna_type. Allowed value: '{TIMEDB_RNA_TYPE}'"
        )


def validate_timedb_expression_type(expression_type: str) -> None:
    if expression_type != TIMEDB_EXPRESSION_TYPE:
        raise DEGPathError(
            f"Invalid TIMEDB expression_type '{expression_type}'. "
            f"Allowed value: '{TIMEDB_EXPRESSION_TYPE}'."
        )


def get_timedb_deg_filename(dataset: str) -> str:
    return f"{dataset}_deg.csv"


def get_timedb_deg_file_path(
    dataset: str,
    rna_type: str,
    expression_type: str,
) -> Path:
    validate_dataset(dataset)
    validate_timedb_rna_type(rna_type)
    validate_timedb_expression_type(expression_type)

    base_dir = get_timedb_dataset_base_dir()

    file_path = (
        base_dir
        / get_timedb_deg_filename(dataset=dataset)
    ).resolve()

    if not str(file_path).startswith(str(base_dir)):
        raise DEGPathError("Invalid TIMEDB DEG file path")

    return file_path


def validate_timedb_deg_file(
    dataset: str,
    rna_type: str,
    expression_type: str,
) -> Path:
    file_path = get_timedb_deg_file_path(
        dataset=dataset,
        rna_type=rna_type,
        expression_type=expression_type,
    )

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"TIMEDB DEG file not found for dataset '{dataset}', "
            f"rna_type '{rna_type}', expression_type '{expression_type}'."
        )

    return file_path


def get_available_timedb_deg_expression_types(
    dataset: str,
    rna_type: str,
) -> list[str]:
    validate_dataset(dataset)
    validate_timedb_rna_type(rna_type)

    file_path = get_timedb_deg_file_path(
        dataset=dataset,
        rna_type=rna_type,
        expression_type=TIMEDB_EXPRESSION_TYPE,
    )

    if file_path.exists() and file_path.is_file():
        return [TIMEDB_EXPRESSION_TYPE]

    return []


def get_tisch2_dataset_base_dir() -> Path:
    return Path(settings.TISCH2_DATASET_BASE_DIR).resolve()


def validate_tisch2_rna_type(rna_type: str) -> None:
    if rna_type != TISCH2_RNA_TYPE:
        raise DEGPathError(
            f"Invalid TISCH2 rna_type. Allowed value: '{TISCH2_RNA_TYPE}'"
        )


def validate_tisch2_expression_type(expression_type: str) -> None:
    if expression_type != TISCH2_EXPRESSION_TYPE:
        raise DEGPathError(
            f"Invalid TISCH2 expression_type '{expression_type}'. "
            f"Allowed value: '{TISCH2_EXPRESSION_TYPE}'."
        )


def get_tisch2_deg_filename(dataset: str) -> str:
    return f"{dataset}_deg.csv"


def get_tisch2_deg_file_path(
    dataset: str,
    rna_type: str,
    expression_type: str,
) -> Path:
    validate_dataset(dataset)
    validate_tisch2_rna_type(rna_type)
    validate_tisch2_expression_type(expression_type)

    base_dir = get_tisch2_dataset_base_dir()

    file_path = (
        base_dir
        / get_tisch2_deg_filename(dataset=dataset)
    ).resolve()

    if not str(file_path).startswith(str(base_dir)):
        raise DEGPathError("Invalid TISCH2 DEG file path")

    return file_path


def validate_tisch2_deg_file(
    dataset: str,
    rna_type: str,
    expression_type: str,
) -> Path:
    file_path = get_tisch2_deg_file_path(
        dataset=dataset,
        rna_type=rna_type,
        expression_type=expression_type,
    )

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"TISCH2 DEG file not found for dataset '{dataset}', "
            f"rna_type '{rna_type}', expression_type '{expression_type}'."
        )

    return file_path


def get_available_tisch2_deg_expression_types(
    dataset: str,
    rna_type: str,
) -> list[str]:
    validate_dataset(dataset)
    validate_tisch2_rna_type(rna_type)

    file_path = get_tisch2_deg_file_path(
        dataset=dataset,
        rna_type=rna_type,
        expression_type=TISCH2_EXPRESSION_TYPE,
    )

    if file_path.exists() and file_path.is_file():
        return [TISCH2_EXPRESSION_TYPE]

    return []


def validate_dataset_deg_file(
    metadata,
    expression_type: str,
) -> Path:
    dataset = metadata.dataset
    rna_type = metadata.gene_bio_type

    if is_tcga_dataset_metadata(metadata):
        return validate_deg_file(
            dataset=dataset,
            rna_type=rna_type,
            expression_type=expression_type,
        )

    if is_cell_dataset_metadata(metadata):
        return validate_tisch2_deg_file(
            dataset=dataset,
            rna_type=rna_type,
            expression_type=expression_type,
        )

    if is_spot_dataset_metadata(metadata):
        raise DEGPathError(
            f"DEG file is not available for spot dataset '{dataset}'."
        )

    return validate_timedb_deg_file(
        dataset=dataset,
        rna_type=rna_type,
        expression_type=expression_type,
    )
