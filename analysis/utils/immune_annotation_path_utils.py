from pathlib import Path

from django.conf import settings


IMMUNE_FILE_PREFIX = "ImmiRImmiR_"
IMMUNE_FILE_SUFFIX = ".csv"


class ImmuneAnnotationPathError(ValueError):
    pass


def get_immune_annotation_raw_dir() -> Path:
    return (Path(settings.IMMUNE_ANNOTATION_BASE_DIR) / "raw").resolve()


def validate_map_info(map_info: str) -> None:
    if not map_info:
        raise ImmuneAnnotationPathError("Missing required parameter: map_info")

    if "/" in map_info or "\\" in map_info or ".." in map_info:
        raise ImmuneAnnotationPathError("Invalid map_info parameter")

    if not map_info.startswith(IMMUNE_FILE_PREFIX):
        raise ImmuneAnnotationPathError(
            f"Invalid map_info. It must start with '{IMMUNE_FILE_PREFIX}'"
        )


def get_label_from_map_info(map_info: str) -> str:
    validate_map_info(map_info)
    return map_info.removeprefix(IMMUNE_FILE_PREFIX)


def get_immune_annotation_file_path(map_info: str) -> Path:
    validate_map_info(map_info)

    raw_dir = get_immune_annotation_raw_dir()
    file_path = (raw_dir / f"{map_info}{IMMUNE_FILE_SUFFIX}").resolve()

    if not str(file_path).startswith(str(raw_dir)):
        raise ImmuneAnnotationPathError("Invalid immune annotation file path")

    return file_path


def validate_immune_annotation_file(map_info: str) -> Path:
    file_path = get_immune_annotation_file_path(map_info)

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Immune annotation file not found for map_info '{map_info}'."
        )

    return file_path
