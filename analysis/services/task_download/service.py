import uuid as uuid_lib
import zipfile
from dataclasses import dataclass
from pathlib import Path

from analysis.services.task_common.registry import (
    find_task_by_uuid,
    TaskNotFoundError,
    MultipleTaskMatchedError,
)
from analysis.services.task_download.registry import (
    get_download_config_by_task_type,
    TaskDownloadConfigNotFoundError,
)
from analysis.services.task_download.resolvers import (
    DownloadableResultFile,
)


class TaskDownloadError(Exception):
    pass


class InvalidTaskUUIDError(TaskDownloadError):
    pass


class TaskNotReadyForDownloadError(TaskDownloadError):
    pass


class TaskResultFileNotFoundError(TaskDownloadError):
    pass


class TaskResultArchiveError(TaskDownloadError):
    pass


@dataclass(frozen=True)
class TaskDownloadResult:
    archive_path: Path
    archive_name: str
    task_type: str
    task_uuid: str


def validate_task_uuid(task_uuid: str) -> None:
    try:
        uuid_lib.UUID(str(task_uuid))
    except ValueError as e:
        raise InvalidTaskUUIDError("Invalid Task UUID format.") from e


def prepare_task_result_download(task_uuid: str) -> TaskDownloadResult:
    """
    Prepare a cached zip archive for downloading workflow task results.

    Rules:
        1. taskUUID must be valid.
        2. Task must exist in the common task registry.
        3. Task status must already be Success.
        4. Download file rules are resolved by task_download.registry.
        5. Required result files must exist.
        6. Optional result files are skipped if missing.
        7. Only files under task output directory are allowed.
        8. Zip archive is generated on demand and cached in output directory.
    """

    task_uuid = str(task_uuid or "").strip()

    if not task_uuid:
        raise InvalidTaskUUIDError("Missing required parameter: taskUUID.")

    validate_task_uuid(task_uuid)

    try:
        task, task_model_config = find_task_by_uuid(task_uuid)
    except TaskNotFoundError as e:
        raise TaskDownloadError(str(e)) from e
    except MultipleTaskMatchedError as e:
        raise TaskDownloadError(str(e)) from e

    if task.status != task.Status.Success:
        raise TaskNotReadyForDownloadError(
            "Task is not completed successfully. "
            f"Current status: {task.get_status_display()}."
        )

    try:
        download_config = get_download_config_by_task_type(
            task_model_config.task_type
        )
    except TaskDownloadConfigNotFoundError as e:
        raise TaskDownloadError(str(e)) from e

    result_files = download_config.result_file_resolver(task)
    archive_name = download_config.archive_name_builder(task)

    output_dir = Path(task.get_output_dir_absolute_path()).resolve()

    downloadable_files = validate_downloadable_result_files(
        result_files=result_files,
        output_dir=output_dir,
    )

    archive_path = get_or_create_cached_result_archive(
        output_dir=output_dir,
        archive_name=archive_name,
        downloadable_files=downloadable_files,
    )

    return TaskDownloadResult(
        archive_path=archive_path,
        archive_name=archive_name,
        task_type=task_model_config.task_type,
        task_uuid=str(task.uuid),
    )


def validate_downloadable_result_files(
    result_files: list[DownloadableResultFile],
    output_dir: Path,
) -> list[DownloadableResultFile]:
    if not output_dir.exists() or not output_dir.is_dir():
        raise TaskResultFileNotFoundError(
            "Task output directory not found."
        )

    if not result_files:
        raise TaskResultFileNotFoundError(
            "No result files are configured for this task."
        )

    downloadable_files = []
    missing_required_files = []

    for result_file in result_files:
        file_path = Path(result_file.path).resolve()

        if not is_path_under_dir(file_path, output_dir):
            raise TaskResultArchiveError(
                f"Invalid result file path: {file_path.name}"
            )

        if file_path.is_symlink():
            raise TaskResultArchiveError(
                "Symbolic link is not allowed in result download: "
                f"{file_path.name}"
            )

        if not file_path.exists() or not file_path.is_file():
            if result_file.required:
                missing_required_files.append(file_path.name)
            continue

        validate_archive_name(result_file.arcname)

        downloadable_files.append(result_file)

    if missing_required_files:
        raise TaskResultFileNotFoundError(
            "Required result file(s) not found: "
            f"{', '.join(missing_required_files)}."
        )

    if not downloadable_files:
        raise TaskResultFileNotFoundError(
            "No downloadable result files found."
        )

    return downloadable_files


def get_or_create_cached_result_archive(
    output_dir: Path,
    archive_name: str,
    downloadable_files: list[DownloadableResultFile],
) -> Path:
    safe_archive_name = sanitize_archive_filename(archive_name)

    archive_path = (output_dir / safe_archive_name).resolve()

    if not is_path_under_dir(archive_path, output_dir):
        raise TaskResultArchiveError(
            "Invalid archive output path."
        )

    if archive_path.exists() and archive_path.is_file():
        return archive_path

    create_cached_result_archive(
        archive_path=archive_path,
        downloadable_files=downloadable_files,
    )

    return archive_path


def create_cached_result_archive(
    archive_path: Path,
    downloadable_files: list[DownloadableResultFile],
) -> None:
    output_dir = archive_path.parent.resolve()

    temp_archive_path = (
        output_dir / f".{archive_path.name}.{uuid_lib.uuid4().hex}.tmp"
    ).resolve()

    if not is_path_under_dir(temp_archive_path, output_dir):
        raise TaskResultArchiveError(
            "Invalid temporary archive output path."
        )

    try:
        with zipfile.ZipFile(
            temp_archive_path,
            mode="w",
            compression=zipfile.ZIP_DEFLATED,
        ) as zip_file:
            for result_file in downloadable_files:
                file_path = Path(result_file.path).resolve()

                if file_path == archive_path:
                    continue

                if file_path == temp_archive_path:
                    continue

                zip_file.write(
                    filename=file_path,
                    arcname=result_file.arcname,
                )

        temp_archive_path.replace(archive_path)

    except Exception as e:
        if temp_archive_path.exists():
            temp_archive_path.unlink()

        raise TaskResultArchiveError(
            f"Failed to create result archive: {str(e)}"
        ) from e


def is_path_under_dir(file_path: Path, base_dir: Path) -> bool:
    try:
        file_path.resolve().relative_to(base_dir.resolve())
        return True
    except ValueError:
        return False


def validate_archive_name(arcname: str) -> None:
    arcname = str(arcname or "").strip()

    if not arcname:
        raise TaskResultArchiveError("Archive entry name cannot be empty.")

    arc_path = Path(arcname)

    if arc_path.is_absolute():
        raise TaskResultArchiveError(
            f"Archive entry name cannot be absolute: {arcname}"
        )

    if ".." in arc_path.parts:
        raise TaskResultArchiveError(
            f"Archive entry name cannot contain '..': {arcname}"
        )


def sanitize_archive_filename(filename: str) -> str:
    filename = str(filename or "").strip()

    if not filename:
        filename = "task_result.zip"

    filename = filename.replace("/", "_").replace("\\", "_").replace("..", "_")

    if not filename.endswith(".zip"):
        filename = f"{filename}.zip"

    return filename
