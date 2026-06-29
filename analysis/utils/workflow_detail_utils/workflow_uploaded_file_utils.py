from pathlib import Path
import uuid as uuid_lib

from analysis.models import PairedCohortTask, HybridReferenceTask
from analysis.utils.paired_cohort_task_utils import (
    PAIRED_COHORT_INPUT_FILENAME_MAP,
    get_paired_cohort_input_file_path,
    PairedCohortTaskPathError,
)
from analysis.utils.hybrid_reference_task_utils import (
    HYBRID_REFERENCE_INPUT_FILENAME_MAP,
    get_hybrid_reference_input_file_path,
    HybridReferenceTaskPathError,
)


class UploadedFileDownloadError(ValueError):
    pass


class InvalidUploadedFileTaskUUIDError(UploadedFileDownloadError):
    pass


class UploadedFileTaskNotFoundError(UploadedFileDownloadError):
    pass


class UploadedFileTypeError(UploadedFileDownloadError):
    pass


class UploadedFilePathError(UploadedFileDownloadError):
    pass


class UploadedFileNotFoundError(UploadedFileDownloadError):
    pass


def validate_uploaded_file_task_uuid(task_uuid: str) -> str:
    task_uuid = str(task_uuid or "").strip()

    if not task_uuid:
        raise InvalidUploadedFileTaskUUIDError(
            "Missing query parameter: taskUUID."
        )

    try:
        uuid_lib.UUID(task_uuid)
    except ValueError as e:
        raise InvalidUploadedFileTaskUUIDError(
            "Invalid Task UUID format."
        ) from e

    return task_uuid


def validate_uploaded_file_type(
    file_type: str,
    filename_map: dict,
) -> str:
    file_type = str(file_type or "").strip()

    if not file_type:
        raise UploadedFileTypeError(
            "Missing query parameter: file_type."
        )

    if file_type not in filename_map:
        raise UploadedFileTypeError(
            "Invalid file_type. Allowed values are: "
            f"{', '.join(filename_map.keys())}."
        )

    return file_type


def get_uploaded_file_response_info(
    *,
    task,
    file_type: str,
    filename_map: dict,
    path_getter,
) -> dict:
    file_type = validate_uploaded_file_type(
        file_type=file_type,
        filename_map=filename_map,
    )

    try:
        file_path = path_getter(
            task=task,
            field_name=file_type,
        )
    except (PairedCohortTaskPathError, HybridReferenceTaskPathError) as e:
        raise UploadedFilePathError(str(e)) from e

    if not file_path.exists() or not file_path.is_file():
        raise UploadedFileNotFoundError(
            "Uploaded file not found: "
            f"{filename_map[file_type]}."
        )

    return {
        "file_path": file_path,
        "filename": file_path.name,
    }


def get_paired_cohort_uploaded_file_response_info(
    *,
    task_uuid: str,
    file_type: str,
) -> dict:
    task_uuid = validate_uploaded_file_task_uuid(task_uuid)

    try:
        task = PairedCohortTask.objects.get(uuid=task_uuid)
    except PairedCohortTask.DoesNotExist as e:
        raise UploadedFileTaskNotFoundError(
            f"PairedCohortTask not found: {task_uuid}."
        ) from e

    return get_uploaded_file_response_info(
        task=task,
        file_type=file_type,
        filename_map=PAIRED_COHORT_INPUT_FILENAME_MAP,
        path_getter=get_paired_cohort_input_file_path,
    )


def get_hybrid_reference_uploaded_file_response_info(
    *,
    task_uuid: str,
    file_type: str,
) -> dict:
    task_uuid = validate_uploaded_file_task_uuid(task_uuid)

    try:
        task = HybridReferenceTask.objects.get(uuid=task_uuid)
    except HybridReferenceTask.DoesNotExist as e:
        raise UploadedFileTaskNotFoundError(
            f"HybridReferenceTask not found: {task_uuid}."
        ) from e

    return get_uploaded_file_response_info(
        task=task,
        file_type=file_type,
        filename_map=HYBRID_REFERENCE_INPUT_FILENAME_MAP,
        path_getter=get_hybrid_reference_input_file_path,
    )
