import traceback

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.services.task_common.registry import TaskNotFoundError, MultipleTaskMatchedError
from analysis.services.task_query.service import (
    query_task_by_uuid,
    InvalidTaskUUIDError,
)
from analysis.services.task_query.status_sync import (
    TaskStatusSyncError,
)


class QueryTaskView(APIView):
    def get(self, request):
        try:
            task_uuid = request.query_params.get("taskUUID", "").strip()

            if not task_uuid:
                return Response(
                    {
                        "success": False,
                        "msg": "Illegal Request! Task UUID is required.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                result = query_task_by_uuid(task_uuid)
            except InvalidTaskUUIDError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )
            except TaskNotFoundError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )
            except MultipleTaskMatchedError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_500_INTERNAL_SERVER_ERROR,
                )
            except TaskStatusSyncError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_500_INTERNAL_SERVER_ERROR,
                )

            return Response(
                result,
                status=status.HTTP_200_OK,
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
