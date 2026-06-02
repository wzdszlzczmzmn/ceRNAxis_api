from pathlib import Path
import os


def get_required_directory(name: str) -> Path:
    value = os.getenv(name)

    if not value:
        raise RuntimeError(
            f"Required environment variable '{name}' is missing."
        )

    path = Path(value).expanduser().resolve()

    if not path.exists():
        raise RuntimeError(
            f"Directory for '{name}' does not exist: {path}"
        )

    if not path.is_dir():
        raise RuntimeError(
            f"Path for '{name}' is not a directory: {path}"
        )

    return path
