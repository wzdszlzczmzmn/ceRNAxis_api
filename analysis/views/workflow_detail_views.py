import traceback
import uuid as uuid_lib
import pandas as pd
import numpy as np

from django.http import FileResponse

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.services.task_common.registry import TaskNotFoundError, MultipleTaskMatchedError
from analysis.services.task_download.registry import TaskDownloadConfigNotFoundError
from analysis.services.task_download.service import prepare_task_result_download, InvalidTaskUUIDError, \
    TaskNotReadyForDownloadError, TaskResultFileNotFoundError, TaskResultArchiveError, TaskDownloadError
from analysis.utils.paired_cohort_task_utils import PairedCohortTaskPathError, PairedCohortTaskInputError, \
    build_paired_cohort_survival_km_data, read_paired_cohort_mrna_gsea_file
from analysis.models import PairedCohortTask


class WorkflowTaskResultDownloadView(APIView):
    """
    Download workflow task result archive.

    Query params:
        taskUUID=<uuid>

    Rules:
        1. Only completed Success tasks can be downloaded.
        2. The result archive is generated on demand.
        3. The generated zip archive is cached under task output directory.
        4. Subsequent downloads reuse the cached archive.
    """

    def get(self, request):
        try:
            task_uuid = str(
                request.query_params.get("taskUUID", "")
            ).strip()

            if not task_uuid:
                return Response(
                    {
                        "success": False,
                        "msg": "Missing required parameter: taskUUID.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                download_result = prepare_task_result_download(
                    task_uuid=task_uuid,
                )

            except InvalidTaskUUIDError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            except TaskNotReadyForDownloadError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            except TaskResultFileNotFoundError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            except TaskResultArchiveError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_500_INTERNAL_SERVER_ERROR,
                )

            except (
                    TaskDownloadError,
                    TaskDownloadConfigNotFoundError,
                    TaskNotFoundError,
                    MultipleTaskMatchedError,
            ) as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            archive_path = download_result.archive_path

            if not archive_path.exists() or not archive_path.is_file():
                return Response(
                    {
                        "success": False,
                        "msg": "Task result archive not found.",
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            return FileResponse(
                open(archive_path, "rb"),
                as_attachment=True,
                filename=download_result.archive_name,
                content_type="application/zip",
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
