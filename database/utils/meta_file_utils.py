from pathlib import Path

from django.conf import settings


class DatasetMetaPathError(ValueError):
    pass


def validate_dataset(dataset: str) -> None:
    if not dataset:
        raise DatasetMetaPathError("Missing required parameter: dataset")

    if "/" in dataset or "\\" in dataset or ".." in dataset:
        raise DatasetMetaPathError("Invalid dataset parameter")


def get_tcga_dataset_meta_file(dataset: str) -> Path:
    validate_dataset(dataset)

    dataset_prefix = "_".join(dataset.split("_")[:2])

    base_dir = Path(settings.TCGA_DATASET_BASE_DIR).resolve()
    file_path = (base_dir / f"{dataset_prefix}_meta.csv").resolve()

    if not str(file_path).startswith(str(base_dir)):
        raise DatasetMetaPathError("Invalid TCGA sample metadata file path")

    return file_path


def get_timedb_dataset_meta_file(dataset: str) -> Path:
    validate_dataset(dataset)

    base_dir = Path(settings.TIMEDB_DATASET_BASE_DIR).resolve()
    file_path = (base_dir / f"{dataset}_meta.csv").resolve()

    if not str(file_path).startswith(str(base_dir)):
        raise DatasetMetaPathError("Invalid TIMEDB sample metadata file path")

    return file_path


def get_dataset_meta_file(dataset: str, expression_mode: str) -> Path:
    if expression_mode == "tcga":
        return get_tcga_dataset_meta_file(dataset)

    if expression_mode == "timedb":
        return get_timedb_dataset_meta_file(dataset)

    raise DatasetMetaPathError(
        f"Invalid expression_mode '{expression_mode}'."
    )
