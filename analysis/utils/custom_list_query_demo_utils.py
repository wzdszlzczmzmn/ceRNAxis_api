import json
from pathlib import Path

from django.conf import settings


CUSTOM_LIST_QUERY_DEMO_DIR_NAME = "custom_list_query"
CUSTOM_LIST_QUERY_DEMO_INPUT_FILENAME = "custom_list_query_demo_input.json"
CUSTOM_LIST_QUERY_DEMO_TASK_TYPE = "CustomListQueryTask"


class CustomListQueryDemoError(ValueError):
    pass


class CustomListQueryDemoPathError(CustomListQueryDemoError):
    pass


class CustomListQueryDemoConfigError(CustomListQueryDemoError):
    pass


def get_custom_list_query_demo_dir() -> Path:
    demo_root = Path(settings.DEMO_INPUT_DATA_HOME).resolve()

    demo_dir = (
        demo_root / CUSTOM_LIST_QUERY_DEMO_DIR_NAME
    ).resolve()

    if not is_path_under_dir(demo_dir, demo_root):
        raise CustomListQueryDemoPathError(
            "Invalid custom list query demo directory."
        )

    return demo_dir


def get_custom_list_query_demo_input_file_path() -> Path:
    demo_dir = get_custom_list_query_demo_dir()

    file_path = (
        demo_dir / CUSTOM_LIST_QUERY_DEMO_INPUT_FILENAME
    ).resolve()

    if not is_path_under_dir(file_path, demo_dir):
        raise CustomListQueryDemoPathError(
            "Invalid custom list query demo input file path."
        )

    return file_path


def load_custom_list_query_demo_input() -> dict:
    file_path = get_custom_list_query_demo_input_file_path()

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Custom list query demo input file not found: {file_path.name}"
        )

    try:
        with file_path.open("r", encoding="utf-8") as f:
            config = json.load(f)
    except json.JSONDecodeError as e:
        raise CustomListQueryDemoConfigError(
            f"Invalid custom list query demo input JSON: {str(e)}"
        ) from e

    validate_custom_list_query_demo_input(config)

    return config


def validate_custom_list_query_demo_input(config: dict) -> None:
    if not isinstance(config, dict):
        raise CustomListQueryDemoConfigError(
            "Custom list query demo input must be a JSON object."
        )

    task_type = str(config.get("task_type", "")).strip()

    if task_type != CUSTOM_LIST_QUERY_DEMO_TASK_TYPE:
        raise CustomListQueryDemoConfigError(
            "Invalid demo task_type. Expected: "
            f"{CUSTOM_LIST_QUERY_DEMO_TASK_TYPE}."
        )

    task_name = str(config.get("task_name", "")).strip()
    map_info = str(config.get("map_info", "")).strip()
    rnas = config.get("rnas")

    if not task_name:
        raise CustomListQueryDemoConfigError(
            "Demo input is missing field: task_name."
        )

    if not map_info:
        raise CustomListQueryDemoConfigError(
            "Demo input is missing field: map_info."
        )

    if not isinstance(rnas, dict):
        raise CustomListQueryDemoConfigError(
            "Demo input field 'rnas' must be an object."
        )


def is_path_under_dir(file_path: Path, base_dir: Path) -> bool:
    try:
        file_path.resolve().relative_to(base_dir.resolve())
        return True
    except ValueError:
        return False
