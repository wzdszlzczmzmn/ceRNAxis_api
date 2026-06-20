import json
import shutil
from pathlib import Path

import pyarrow.parquet as pq
from django.conf import settings

from analysis.utils.paired_cohort_task_utils import (
    PAIRED_COHORT_ALLOWED_FILE_FIELDS,
    get_paired_cohort_input_file_path,
    PairedCohortTaskInputError,
)


PAIRED_COHORT_DEMO_WORKFLOW_TYPE = "paired_cohort"

PAIRED_COHORT_DEMO_MANIFEST_FILENAME = "demo_manifest.json"

PAIRED_COHORT_DEMO_DATA_ARCHIVE_NAME = "paired_cohort_demo_data.zip"

PAIRED_COHORT_DEMO_RNA_TYPE_TO_FILE_KEY = {
    "mRNA": "mrna_file",
    "miRNA": "mirna_file",
    "lncRNA": "lncrna_file",
}

PAIRED_COHORT_DEMO_INPUT_FILENAME = "paired_cohort_demo_input.json"

PAIRED_COHORT_DEMO_TASK_TYPE = "PairedCohortTask"

PAIRED_COHORT_DEMO_ALLOWED_DEG_METHODS = [
    "limma",
    "deseq2",
]

PAIRED_COHORT_DEMO_CUTOFF_RNA_TYPES = [
    "mRNA",
    "miRNA",
    "lncRNA",
]

PAIRED_COHORT_DEMO_VALID_RNA_TYPES = list(
    PAIRED_COHORT_DEMO_RNA_TYPE_TO_FILE_KEY.keys()
)


class PairedCohortDemoError(ValueError):
    pass


class PairedCohortDemoPathError(PairedCohortDemoError):
    pass


class PairedCohortDemoManifestError(PairedCohortDemoError):
    pass


class PairedCohortDemoConfigError(PairedCohortDemoError):
    pass


def get_paired_cohort_demo_dir() -> Path:
    demo_root = Path(settings.DEMO_INPUT_DATA_HOME).resolve()
    demo_dir = (demo_root / PAIRED_COHORT_DEMO_WORKFLOW_TYPE).resolve()

    if not is_path_under_dir(demo_dir, demo_root):
        raise PairedCohortDemoPathError(
            "Invalid paired cohort demo directory."
        )

    return demo_dir


def get_paired_cohort_demo_manifest_path() -> Path:
    demo_dir = get_paired_cohort_demo_dir()

    manifest_path = (
        demo_dir / PAIRED_COHORT_DEMO_MANIFEST_FILENAME
    ).resolve()

    if not is_path_under_dir(manifest_path, demo_dir):
        raise PairedCohortDemoPathError(
            "Invalid paired cohort demo manifest path."
        )

    return manifest_path


def load_paired_cohort_demo_manifest() -> dict:
    manifest_path = get_paired_cohort_demo_manifest_path()

    if not manifest_path.exists() or not manifest_path.is_file():
        raise FileNotFoundError(
            f"Paired cohort demo manifest not found: {manifest_path.name}"
        )

    try:
        with manifest_path.open("r", encoding="utf-8") as f:
            manifest = json.load(f)
    except json.JSONDecodeError as e:
        raise PairedCohortDemoManifestError(
            f"Invalid paired cohort demo manifest JSON: {str(e)}"
        ) from e

    validate_paired_cohort_demo_manifest(manifest)

    return manifest


def validate_paired_cohort_demo_manifest(manifest: dict) -> None:
    if not isinstance(manifest, dict):
        raise PairedCohortDemoManifestError(
            "Paired cohort demo manifest must be a JSON object."
        )

    csv_files = manifest.get("csv_files")
    parquet_files = manifest.get("parquet_files")

    if not isinstance(csv_files, dict):
        raise PairedCohortDemoManifestError(
            "Manifest field 'csv_files' must be an object."
        )

    if not isinstance(parquet_files, dict):
        raise PairedCohortDemoManifestError(
            "Manifest field 'parquet_files' must be an object."
        )

    required_csv_keys = [
        "mrna_file",
        "mirna_file",
        "lncrna_file",
        "meta_file",
    ]

    required_parquet_keys = [
        "mrna_file",
        "mirna_file",
        "lncrna_file",
    ]

    missing_csv_keys = [
        key for key in required_csv_keys
        if key not in csv_files
    ]

    missing_parquet_keys = [
        key for key in required_parquet_keys
        if key not in parquet_files
    ]

    if missing_csv_keys:
        raise PairedCohortDemoManifestError(
            "Manifest field 'csv_files' is missing key(s): "
            f"{', '.join(missing_csv_keys)}."
        )

    if missing_parquet_keys:
        raise PairedCohortDemoManifestError(
            "Manifest field 'parquet_files' is missing key(s): "
            f"{', '.join(missing_parquet_keys)}."
        )


def get_demo_file_path_from_manifest(
    manifest: dict,
    file_group: str,
    file_key: str,
) -> Path:
    if file_group not in ["csv_files", "parquet_files"]:
        raise PairedCohortDemoPathError(
            f"Invalid demo file group: {file_group}."
        )

    files = manifest.get(file_group)

    if not isinstance(files, dict):
        raise PairedCohortDemoManifestError(
            f"Manifest field '{file_group}' must be an object."
        )

    filename = str(files.get(file_key, "")).strip()

    if not filename:
        raise PairedCohortDemoManifestError(
            f"Manifest field '{file_group}.{file_key}' is missing."
        )

    validate_demo_filename(filename)

    demo_dir = get_paired_cohort_demo_dir()
    file_path = (demo_dir / filename).resolve()

    if not is_path_under_dir(file_path, demo_dir):
        raise PairedCohortDemoPathError(
            f"Invalid demo file path: {filename}."
        )

    return file_path


def get_paired_cohort_demo_meta_file_path() -> Path:
    manifest = load_paired_cohort_demo_manifest()

    return get_demo_file_path_from_manifest(
        manifest=manifest,
        file_group="csv_files",
        file_key="meta_file",
    )


def get_paired_cohort_demo_expression_file_path(
    rna_type: str,
) -> Path:
    rna_type = str(rna_type or "").strip()

    if rna_type not in PAIRED_COHORT_DEMO_RNA_TYPE_TO_FILE_KEY:
        raise PairedCohortDemoManifestError(
            "Invalid rna_type. Allowed values are: "
            f"{', '.join(PAIRED_COHORT_DEMO_VALID_RNA_TYPES)}."
        )

    file_key = PAIRED_COHORT_DEMO_RNA_TYPE_TO_FILE_KEY[rna_type]

    manifest = load_paired_cohort_demo_manifest()

    return get_demo_file_path_from_manifest(
        manifest=manifest,
        file_group="parquet_files",
        file_key=file_key,
    )


def validate_demo_file_exists(file_path: Path, file_label: str) -> Path:
    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"{file_label} not found: {file_path.name}"
        )

    return file_path


def read_parquet_columns(file_path: Path) -> list[str]:
    try:
        schema = pq.read_schema(file_path)
        return schema.names
    except Exception as e:
        raise PairedCohortDemoError(
            f"Failed to read parquet schema: {str(e)}"
        ) from e


def validate_demo_filename(filename: str) -> None:
    filename = str(filename or "").strip()

    if not filename:
        raise PairedCohortDemoPathError(
            "Demo filename cannot be empty."
        )

    path = Path(filename)

    if path.is_absolute():
        raise PairedCohortDemoPathError(
            f"Demo filename cannot be absolute: {filename}"
        )

    if ".." in path.parts:
        raise PairedCohortDemoPathError(
            f"Demo filename cannot contain '..': {filename}"
        )

    if "/" in filename or "\\" in filename:
        raise PairedCohortDemoPathError(
            f"Demo filename cannot contain path separators: {filename}"
        )


def is_path_under_dir(file_path: Path, base_dir: Path) -> bool:
    try:
        file_path.resolve().relative_to(base_dir.resolve())
        return True
    except ValueError:
        return False


def get_paired_cohort_demo_data_archive_path() -> Path:
    demo_dir = get_paired_cohort_demo_dir()

    archive_path = (
        demo_dir / PAIRED_COHORT_DEMO_DATA_ARCHIVE_NAME
    ).resolve()

    if not is_path_under_dir(archive_path, demo_dir):
        raise PairedCohortDemoPathError(
            "Invalid paired cohort demo data archive path."
        )

    return archive_path


def validate_paired_cohort_demo_data_archive() -> Path:
    archive_path = get_paired_cohort_demo_data_archive_path()

    if not archive_path.exists() or not archive_path.is_file():
        raise FileNotFoundError(
            f"Paired cohort demo data archive not found: {archive_path.name}"
        )

    return archive_path


def get_paired_cohort_demo_input_config_path() -> Path:
    demo_dir = get_paired_cohort_demo_dir()

    config_path = (
        demo_dir / PAIRED_COHORT_DEMO_INPUT_FILENAME
    ).resolve()

    if not is_path_under_dir(config_path, demo_dir):
        raise PairedCohortDemoPathError(
            "Invalid paired cohort demo input config path."
        )

    return config_path


def load_paired_cohort_demo_input() -> dict:
    config_path = get_paired_cohort_demo_input_config_path()

    if not config_path.exists() or not config_path.is_file():
        raise FileNotFoundError(
            f"Paired cohort demo input config not found: {config_path.name}"
        )

    try:
        with config_path.open("r", encoding="utf-8") as f:
            config = json.load(f)
    except json.JSONDecodeError as e:
        raise PairedCohortDemoConfigError(
            f"Invalid paired cohort demo input JSON: {str(e)}"
        ) from e

    validate_paired_cohort_demo_input(config)

    return config


def validate_paired_cohort_demo_input(config: dict) -> None:
    if not isinstance(config, dict):
        raise PairedCohortDemoConfigError(
            "Paired cohort demo input must be a JSON object."
        )

    task_type = str(config.get("task_type", "")).strip()

    if task_type != PAIRED_COHORT_DEMO_TASK_TYPE:
        raise PairedCohortDemoConfigError(
            f"Invalid demo task_type. Expected: {PAIRED_COHORT_DEMO_TASK_TYPE}."
        )

    task_name = str(config.get("task_name", "")).strip()
    map_info = str(config.get("map_info", "")).strip()
    deg_method = str(config.get("deg_method", "")).strip()

    if not task_name:
        raise PairedCohortDemoConfigError(
            "Demo input is missing field: task_name."
        )

    if not map_info:
        raise PairedCohortDemoConfigError(
            "Demo input is missing field: map_info."
        )

    if deg_method not in PAIRED_COHORT_DEMO_ALLOWED_DEG_METHODS:
        raise PairedCohortDemoConfigError(
            "Invalid demo deg_method. Allowed values are: "
            f"{', '.join(PAIRED_COHORT_DEMO_ALLOWED_DEG_METHODS)}."
        )

    files = config.get("files")

    if not isinstance(files, dict):
        raise PairedCohortDemoConfigError(
            "Demo input field 'files' must be an object."
        )

    missing_files = [
        field_name
        for field_name in PAIRED_COHORT_ALLOWED_FILE_FIELDS
        if field_name not in files
    ]

    if missing_files:
        raise PairedCohortDemoConfigError(
            "Demo input field 'files' is missing key(s): "
            f"{', '.join(missing_files)}."
        )

    cutoffs = config.get("cutoffs")

    if not isinstance(cutoffs, dict):
        raise PairedCohortDemoConfigError(
            "Demo input field 'cutoffs' must be an object."
        )

    for rna_type in PAIRED_COHORT_DEMO_CUTOFF_RNA_TYPES:
        cutoff_config = cutoffs.get(rna_type)

        if not isinstance(cutoff_config, dict):
            raise PairedCohortDemoConfigError(
                f"Demo input field 'cutoffs.{rna_type}' must be an object."
            )

        for cutoff_name in ["logfc_cutoff", "padj_cutoff"]:
            if cutoff_name not in cutoff_config:
                raise PairedCohortDemoConfigError(
                    f"Demo input field 'cutoffs.{rna_type}.{cutoff_name}' is missing."
                )

            parse_demo_cutoff(
                config=config,
                rna_type=rna_type,
                cutoff_name=cutoff_name,
            )


def get_paired_cohort_demo_source_file_path(
    config: dict,
    field_name: str,
) -> Path:
    if field_name not in PAIRED_COHORT_ALLOWED_FILE_FIELDS:
        raise PairedCohortDemoConfigError(
            f"Invalid paired cohort demo file field: {field_name}."
        )

    files = config.get("files")

    if not isinstance(files, dict):
        raise PairedCohortDemoConfigError(
            "Demo input field 'files' must be an object."
        )

    filename = str(files.get(field_name, "")).strip()

    if not filename:
        raise PairedCohortDemoConfigError(
            f"Demo input field 'files.{field_name}' is missing."
        )

    validate_demo_filename(filename)

    demo_dir = get_paired_cohort_demo_dir()
    file_path = (demo_dir / filename).resolve()

    if not is_path_under_dir(file_path, demo_dir):
        raise PairedCohortDemoPathError(
            f"Invalid paired cohort demo source file path: {filename}."
        )

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Paired cohort demo source file not found: {filename}"
        )

    return file_path


def copy_paired_cohort_demo_input_files_to_task(
    task,
    config: dict,
) -> dict:
    """
    Copy paired cohort demo CSV files into task workspace/input.

    Source files are defined in paired_cohort_demo_input.json.

    Destination filenames are controlled by PAIRED_COHORT_INPUT_FILENAME_MAP
    through get_paired_cohort_input_file_path(), for example:

        mrna_file   -> input/mrna.csv
        mirna_file  -> input/mirna.csv
        lncrna_file -> input/lncrna.csv
        meta_file   -> input/meta.csv
    """

    saved_files = {}

    for field_name in PAIRED_COHORT_ALLOWED_FILE_FIELDS:
        source_path = get_paired_cohort_demo_source_file_path(
            config=config,
            field_name=field_name,
        )

        dest_path = get_paired_cohort_input_file_path(
            task=task,
            field_name=field_name,
        )

        dest_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            shutil.copyfile(source_path, dest_path)
        except OSError as e:
            raise PairedCohortTaskInputError(
                f"Failed to copy demo input file '{source_path.name}': {str(e)}"
            ) from e

        saved_files[field_name] = dest_path.name

    return saved_files


def parse_demo_cutoff(
    config: dict,
    rna_type: str,
    cutoff_name: str,
) -> float:
    try:
        raw_value = config["cutoffs"][rna_type][cutoff_name]
        value = float(raw_value)
    except (KeyError, TypeError, ValueError) as e:
        raise PairedCohortDemoConfigError(
            f"Invalid cutoff value: cutoffs.{rna_type}.{cutoff_name}."
        ) from e

    if cutoff_name == "logfc_cutoff" and value < 0:
        raise PairedCohortDemoConfigError(
            f"cutoffs.{rna_type}.{cutoff_name} must be greater than or equal to 0."
        )

    if cutoff_name == "padj_cutoff" and (value <= 0 or value > 1):
        raise PairedCohortDemoConfigError(
            f"cutoffs.{rna_type}.{cutoff_name} must be in the range (0, 1]."
        )

    return value


def get_paired_cohort_demo_cutoff_fields(config: dict) -> dict:
    """
    Return cutoff fields matching PairedCohortTask model field names.
    """

    return {
        "logfc_cutoff_mrna": parse_demo_cutoff(
            config=config,
            rna_type="mRNA",
            cutoff_name="logfc_cutoff",
        ),
        "padj_cutoff_mrna": parse_demo_cutoff(
            config=config,
            rna_type="mRNA",
            cutoff_name="padj_cutoff",
        ),
        "logfc_cutoff_mirna": parse_demo_cutoff(
            config=config,
            rna_type="miRNA",
            cutoff_name="logfc_cutoff",
        ),
        "padj_cutoff_mirna": parse_demo_cutoff(
            config=config,
            rna_type="miRNA",
            cutoff_name="padj_cutoff",
        ),
        "logfc_cutoff_lncrna": parse_demo_cutoff(
            config=config,
            rna_type="lncRNA",
            cutoff_name="logfc_cutoff",
        ),
        "padj_cutoff_lncrna": parse_demo_cutoff(
            config=config,
            rna_type="lncRNA",
            cutoff_name="padj_cutoff",
        ),
    }
