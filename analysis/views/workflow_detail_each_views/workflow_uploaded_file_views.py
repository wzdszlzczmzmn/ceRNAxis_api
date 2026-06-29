import traceback

from django.http import FileResponse

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.utils.workflow_detail_utils.workflow_uploaded_file_utils import (
    UploadedFileDownloadError,
    InvalidUploadedFileTaskUUIDError,
    UploadedFileTaskNotFoundError,
    UploadedFileTypeError,
    UploadedFilePathError,
    UploadedFileNotFoundError,
    get_paired_cohort_uploaded_file_response_info,
    get_hybrid_reference_uploaded_file_response_info,
)


class BaseUploadedFileDownloadView(APIView):
    """
    Base view for downloading workflow uploaded input files.

    Subclasses should define:
        get_file_response_info(task_uuid, file_type)
    """

    def get(self, request):
        try:
            task_uuid = str(
                request.query_params.get("taskUUID", "")
            ).strip()

            file_type = str(
                request.query_params.get("file_type", "")
            ).strip()

            file_info = self.get_file_response_info(
                task_uuid=task_uuid,
                file_type=file_type,
            )

            file_path = file_info["file_path"]
            filename = file_info["filename"]

            return FileResponse(
                open(file_path, "rb"),
                as_attachment=True,
                filename=filename,
                content_type="text/csv",
            )

        except InvalidUploadedFileTaskUUIDError as e:
            return Response(
                {
                    "success": False,
                    "msg": str(e),
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        except UploadedFileTypeError as e:
            return Response(
                {
                    "success": False,
                    "msg": str(e),
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        except UploadedFilePathError as e:
            return Response(
                {
                    "success": False,
                    "msg": str(e),
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        except UploadedFileTaskNotFoundError as e:
            return Response(
                {
                    "success": False,
                    "msg": str(e),
                },
                status=status.HTTP_404_NOT_FOUND,
            )

        except UploadedFileNotFoundError as e:
            return Response(
                {
                    "success": False,
                    "msg": str(e),
                },
                status=status.HTTP_404_NOT_FOUND,
            )

        except UploadedFileDownloadError as e:
            return Response(
                {
                    "success": False,
                    "msg": str(e),
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {
                    "success": False,
                    "msg": f"Server error: {str(e)}",
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

    def get_file_response_info(self, *, task_uuid: str, file_type: str) -> dict:
        raise NotImplementedError


class PairedCohortUploadedFileDownloadView(BaseUploadedFileDownloadView):
    """
    Download uploaded input file for PairedCohortTask.

    Query params:
        taskUUID: PairedCohortTask UUID
        file_type: one of mrna_file, mirna_file, lncrna_file, circrna_file, meta_file
    """

    def get_file_response_info(self, *, task_uuid: str, file_type: str) -> dict:
        return get_paired_cohort_uploaded_file_response_info(
            task_uuid=task_uuid,
            file_type=file_type,
        )


class HybridReferenceUploadedFileDownloadView(BaseUploadedFileDownloadView):
    """
    Download uploaded input file for HybridReferenceTask.

    Query params:
        taskUUID: HybridReferenceTask UUID
        file_type: one of mrna_file, meta_file
    """

    def get_file_response_info(self, *, task_uuid: str, file_type: str) -> dict:
        return get_hybrid_reference_uploaded_file_response_info(
            task_uuid=task_uuid,
            file_type=file_type,
        )
