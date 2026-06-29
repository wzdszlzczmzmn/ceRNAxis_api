import json
import shutil
from pathlib import Path
import pyarrow.parquet as pq

from django.conf import settings

from analysis.utils.hybrid_reference_task_utils import (
    HYBRID_REFERENCE_ALLOWED_FILE_FIELDS,
    HYBRID_REFERENCE_VALID_DEG_METHODS,
    HYBRID_REFERENCE_VALID_LNCRNA_TYPES,
    HYBRID_REFERENCE_VALID_TCGA_TYPES,
    get_hybrid_reference_input_file_path,
)


HYBRID_REFERENCE_DEMO_WORKFLOW_TYPE = "hybrid_reference"

HYBRID_REFERENCE_DEMO_MANIFEST_FILENAME = "demo_manifest.json"

HYBRID_REFERENCE_DEMO_DATA_ARCHIVE_NAME = "hybrid_reference_demo_data.zip"

HYBRID_REFERENCE_DEMO_INPUT_FILENAME = "hybrid_reference_demo_input.json"

HYBRID_REFERENCE_DEMO_TASK_TYPE = "HybridReferenceTask"

HYBRID_REFERENCE_DEMO_CUTOFF_RNA_TYPES = [
    "mRNA",
]

HYBRID_REFERENCE_DEMO_RNA_TYPE_TO_FILE_KEY = {
    "mRNA": "mrna_file",
}

HYBRID_REFERENCE_DEMO_VALID_RNA_TYPES = list(
    HYBRID_REFERENCE_DEMO_RNA_TYPE_TO_FILE_KEY.keys()
)


class HybridReferenceDemoError(ValueError):
    pass


class HybridReferenceDemoPathError(HybridReferenceDemoError):
    pass


class HybridReferenceDemoManifestError(HybridReferenceDemoError):
    pass


class HybridReferenceDemoConfigError(HybridReferenceDemoError):
    pass


def is_path_under_dir(file_path: Path, base_dir: Path) -> bool:
    try:
        file_path.resolve().relative_to(base_dir.resolve())
        return True
    except ValueError:
        return False


def validate_demo_filename(filename: str) -> None:
    filename = str(filename or "").strip()

    if not filename:
        raise HybridReferenceDemoPathError(
            "Demo filename cannot be empty."
        )

    path = Path(filename)

    if path.is_absolute():
        raise HybridReferenceDemoPathError(
            f"Demo filename cannot be absolute: {filename}"
        )

    if ".." in path.parts:
        raise HybridReferenceDemoPathError(
            f"Demo filename cannot contain '..': {filename}"
        )

    if "/" in filename or "\\" in filename:
        raise HybridReferenceDemoPathError(
            f"Demo filename cannot contain path separators: {filename}"
        )


def get_hybrid_reference_demo_dir() -> Path:
    demo_root = Path(settings.DEMO_INPUT_DATA_HOME).resolve()
    demo_dir = (demo_root / HYBRID_REFERENCE_DEMO_WORKFLOW_TYPE).resolve()

    if not is_path_under_dir(demo_dir, demo_root):
        raise HybridReferenceDemoPathError(
            "Invalid hybrid reference demo directory."
        )

    return demo_dir


def get_hybrid_reference_demo_manifest_path() -> Path:
    demo_dir = get_hybrid_reference_demo_dir()

    manifest_path = (
        demo_dir / HYBRID_REFERENCE_DEMO_MANIFEST_FILENAME
    ).resolve()

    if not is_path_under_dir(manifest_path, demo_dir):
        raise HybridReferenceDemoPathError(
            "Invalid hybrid reference demo manifest path."
        )

    return manifest_path


def validate_hybrid_reference_demo_manifest(manifest: dict) -> None:
    if not isinstance(manifest, dict):
        raise HybridReferenceDemoManifestError(
            "Hybrid reference demo manifest must be a JSON object."
        )

    task_type = str(manifest.get("task_type", "")).strip()

    if task_type != HYBRID_REFERENCE_DEMO_TASK_TYPE:
        raise HybridReferenceDemoManifestError(
            f"Invalid demo task_type. Expected: {HYBRID_REFERENCE_DEMO_TASK_TYPE}."
        )

    csv_files = manifest.get("csv_files")
    parquet_files = manifest.get("parquet_files")

    if not isinstance(csv_files, dict):
        raise HybridReferenceDemoManifestError(
            "Manifest field 'csv_files' must be an object."
        )

    if not isinstance(parquet_files, dict):
        raise HybridReferenceDemoManifestError(
            "Manifest field 'parquet_files' must be an object."
        )

    required_csv_keys = [
        "mrna_file",
        "meta_file",
    ]

    required_parquet_keys = [
        "mrna_file",
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
        raise HybridReferenceDemoManifestError(
            "Manifest field 'csv_files' is missing key(s): "
            f"{', '.join(missing_csv_keys)}."
        )

    if missing_parquet_keys:
        raise HybridReferenceDemoManifestError(
            "Manifest field 'parquet_files' is missing key(s): "
            f"{', '.join(missing_parquet_keys)}."
        )


def load_hybrid_reference_demo_manifest() -> dict:
    manifest_path = get_hybrid_reference_demo_manifest_path()

    if not manifest_path.exists() or not manifest_path.is_file():
        raise FileNotFoundError(
            f"Hybrid reference demo manifest not found: {manifest_path.name}"
        )

    try:
        with manifest_path.open("r", encoding="utf-8") as f:
            manifest = json.load(f)
    except json.JSONDecodeError as e:
        raise HybridReferenceDemoManifestError(
            f"Invalid hybrid reference demo manifest JSON: {str(e)}"
        ) from e

    validate_hybrid_reference_demo_manifest(manifest)

    return manifest


def get_hybrid_reference_demo_input_config_path() -> Path:
    demo_dir = get_hybrid_reference_demo_dir()

    config_path = (
        demo_dir / HYBRID_REFERENCE_DEMO_INPUT_FILENAME
    ).resolve()

    if not is_path_under_dir(config_path, demo_dir):
        raise HybridReferenceDemoPathError(
            "Invalid hybrid reference demo input config path."
        )

    return config_path


def load_hybrid_reference_demo_input() -> dict:
    config_path = get_hybrid_reference_demo_input_config_path()

    if not config_path.exists() or not config_path.is_file():
        raise FileNotFoundError(
            f"Hybrid reference demo input config not found: {config_path.name}"
        )

    try:
        with config_path.open("r", encoding="utf-8") as f:
            config = json.load(f)
    except json.JSONDecodeError as e:
        raise HybridReferenceDemoConfigError(
            f"Invalid hybrid reference demo input JSON: {str(e)}"
        ) from e

    validate_hybrid_reference_demo_input(config)

    return config


def validate_hybrid_reference_demo_input(config: dict) -> None:
    if not isinstance(config, dict):
        raise HybridReferenceDemoConfigError(
            "Hybrid reference demo input must be a JSON object."
        )

    task_type = str(config.get("task_type", "")).strip()

    if task_type != HYBRID_REFERENCE_DEMO_TASK_TYPE:
        raise HybridReferenceDemoConfigError(
            f"Invalid demo task_type. Expected: {HYBRID_REFERENCE_DEMO_TASK_TYPE}."
        )

    task_name = str(config.get("task_name", "")).strip()
    map_info = str(config.get("map_info", "")).strip()
    tcga_type = str(config.get("tcga_type", "")).strip()
    lncrna_type = str(config.get("lncrna_type", "")).strip()
    deg_method = str(config.get("deg_method", "")).strip()
    use_padj = config.get("use_padj", True)

    if not task_name:
        raise HybridReferenceDemoConfigError(
            "Demo input is missing field: task_name."
        )

    if not map_info:
        raise HybridReferenceDemoConfigError(
            "Demo input is missing field: map_info."
        )

    if tcga_type not in HYBRID_REFERENCE_VALID_TCGA_TYPES:
        raise HybridReferenceDemoConfigError(
            "Invalid demo tcga_type. Allowed values are: "
            f"{', '.join(HYBRID_REFERENCE_VALID_TCGA_TYPES)}."
        )

    if lncrna_type not in HYBRID_REFERENCE_VALID_LNCRNA_TYPES:
        raise HybridReferenceDemoConfigError(
            "Invalid demo lncrna_type. Allowed values are: "
            f"{', '.join(HYBRID_REFERENCE_VALID_LNCRNA_TYPES)}."
        )

    if deg_method not in HYBRID_REFERENCE_VALID_DEG_METHODS:
        raise HybridReferenceDemoConfigError(
            "Invalid demo deg_method. Allowed values are: "
            f"{', '.join(HYBRID_REFERENCE_VALID_DEG_METHODS)}."
        )

    if not isinstance(use_padj, bool):
        raise HybridReferenceDemoConfigError(
            "Demo input field 'use_padj' must be a boolean."
        )

    files = config.get("files")

    if not isinstance(files, dict):
        raise HybridReferenceDemoConfigError(
            "Demo input field 'files' must be an object."
        )

    missing_files = [
        field_name
        for field_name in HYBRID_REFERENCE_ALLOWED_FILE_FIELDS
        if field_name not in files
    ]

    if missing_files:
        raise HybridReferenceDemoConfigError(
            "Demo input field 'files' is missing key(s): "
            f"{', '.join(missing_files)}."
        )

    cutoffs = config.get("cutoffs")

    if not isinstance(cutoffs, dict):
        raise HybridReferenceDemoConfigError(
            "Demo input field 'cutoffs' must be an object."
        )

    for rna_type in HYBRID_REFERENCE_DEMO_CUTOFF_RNA_TYPES:
        cutoff_config = cutoffs.get(rna_type)

        if not isinstance(cutoff_config, dict):
            raise HybridReferenceDemoConfigError(
                f"Demo input field 'cutoffs.{rna_type}' must be an object."
            )

        for cutoff_name in ["logfc_cutoff", "padj_cutoff"]:
            if cutoff_name not in cutoff_config:
                raise HybridReferenceDemoConfigError(
                    f"Demo input field 'cutoffs.{rna_type}.{cutoff_name}' is missing."
                )

            parse_demo_cutoff(
                config=config,
                rna_type=rna_type,
                cutoff_name=cutoff_name,
            )


def get_hybrid_reference_demo_source_file_path(
    config: dict,
    field_name: str,
) -> Path:
    if field_name not in HYBRID_REFERENCE_ALLOWED_FILE_FIELDS:
        raise HybridReferenceDemoConfigError(
            f"Invalid hybrid reference demo file field: {field_name}."
        )

    files = config.get("files")

    if not isinstance(files, dict):
        raise HybridReferenceDemoConfigError(
            "Demo input field 'files' must be an object."
        )

    filename = str(files.get(field_name, "")).strip()

    if not filename:
        raise HybridReferenceDemoConfigError(
            f"Demo input field 'files.{field_name}' is missing."
        )

    validate_demo_filename(filename)

    demo_dir = get_hybrid_reference_demo_dir()
    file_path = (demo_dir / filename).resolve()

    if not is_path_under_dir(file_path, demo_dir):
        raise HybridReferenceDemoPathError(
            f"Invalid hybrid reference demo source file path: {filename}."
        )

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Hybrid reference demo source file not found: {filename}"
        )

    return file_path


def copy_hybrid_reference_demo_input_files_to_task(
    *,
    task,
    config: dict,
) -> dict:
    saved_files = {}

    for field_name in HYBRID_REFERENCE_ALLOWED_FILE_FIELDS:
        source_file_path = get_hybrid_reference_demo_source_file_path(
            config=config,
            field_name=field_name,
        )

        target_file_path = get_hybrid_reference_input_file_path(
            task=task,
            field_name=field_name,
        )

        target_file_path.parent.mkdir(parents=True, exist_ok=True)

        shutil.copyfile(source_file_path, target_file_path)

        saved_files[field_name] = target_file_path.name

    return saved_files


def parse_demo_cutoff(
    *,
    config: dict,
    rna_type: str,
    cutoff_name: str,
) -> float:
    try:
        value = float(config["cutoffs"][rna_type][cutoff_name])
    except (KeyError, TypeError, ValueError) as e:
        raise HybridReferenceDemoConfigError(
            f"Invalid cutoff value: cutoffs.{rna_type}.{cutoff_name}."
        ) from e

    if cutoff_name == "logfc_cutoff" and value < 0:
        raise HybridReferenceDemoConfigError(
            f"cutoffs.{rna_type}.{cutoff_name} must be greater than or equal to 0."
        )

    if cutoff_name == "padj_cutoff" and (value <= 0 or value > 1):
        raise HybridReferenceDemoConfigError(
            f"cutoffs.{rna_type}.{cutoff_name} must be in the range (0, 1]."
        )

    return value


def get_hybrid_reference_demo_cutoff_fields(config: dict) -> dict:
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
    }


def get_hybrid_reference_demo_data_archive_path() -> Path:
    demo_dir = get_hybrid_reference_demo_dir()

    archive_path = (
        demo_dir / HYBRID_REFERENCE_DEMO_DATA_ARCHIVE_NAME
    ).resolve()

    if not is_path_under_dir(archive_path, demo_dir):
        raise HybridReferenceDemoPathError(
            "Invalid hybrid reference demo data archive path."
        )

    return archive_path


def validate_hybrid_reference_demo_data_archive() -> Path:
    archive_path = get_hybrid_reference_demo_data_archive_path()

    if not archive_path.exists() or not archive_path.is_file():
        raise FileNotFoundError(
            f"Hybrid reference demo data archive not found: {archive_path.name}"
        )

    return archive_path


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
        raise HybridReferenceDemoError(
            f"Failed to read parquet schema: {str(e)}"
        ) from e


def get_demo_file_path_from_manifest(
    manifest: dict,
    file_group: str,
    file_key: str,
) -> Path:
    if file_group not in ["csv_files", "parquet_files"]:
        raise HybridReferenceDemoPathError(
            f"Invalid demo file group: {file_group}."
        )

    files = manifest.get(file_group)

    if not isinstance(files, dict):
        raise HybridReferenceDemoManifestError(
            f"Manifest field '{file_group}' must be an object."
        )

    filename = str(files.get(file_key, "")).strip()

    if not filename:
        raise HybridReferenceDemoManifestError(
            f"Manifest field '{file_group}.{file_key}' is missing."
        )

    validate_demo_filename(filename)

    demo_dir = get_hybrid_reference_demo_dir()
    file_path = (demo_dir / filename).resolve()

    if not is_path_under_dir(file_path, demo_dir):
        raise HybridReferenceDemoPathError(
            f"Invalid hybrid reference demo file path: {filename}."
        )

    return file_path


def get_hybrid_reference_demo_meta_file_path() -> Path:
    manifest = load_hybrid_reference_demo_manifest()

    return get_demo_file_path_from_manifest(
        manifest=manifest,
        file_group="csv_files",
        file_key="meta_file",
    )


def get_hybrid_reference_demo_expression_file_path(
    rna_type: str,
) -> Path:
    rna_type = str(rna_type or "").strip()

    if rna_type not in HYBRID_REFERENCE_DEMO_RNA_TYPE_TO_FILE_KEY:
        raise HybridReferenceDemoManifestError(
            "Invalid rna_type. Allowed values are: "
            f"{', '.join(HYBRID_REFERENCE_DEMO_VALID_RNA_TYPES)}."
        )

    file_key = HYBRID_REFERENCE_DEMO_RNA_TYPE_TO_FILE_KEY[rna_type]

    manifest = load_hybrid_reference_demo_manifest()

    return get_demo_file_path_from_manifest(
        manifest=manifest,
        file_group="parquet_files",
        file_key=file_key,
    )
