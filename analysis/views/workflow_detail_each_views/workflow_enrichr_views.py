import traceback

from rest_framework import status
from rest_framework.response import Response
from rest_framework.views import APIView

from analysis.models import CustomListQueryTask

from analysis.utils.workflow_detail_utils.workflow_network_view_utils import (
    WorkflowNetworkViewError,
    get_required_task_uuid,
    get_task_or_raise,
    require_success_task,
)

from analysis.utils.workflow_detail_utils.workflow_enrichr_utils import (
    CustomListQueryEnrichrInputError,
    CustomListQueryEnrichrPathError,
    normalize_enrichr_direction,
    read_custom_list_query_enrichr_file,
    build_custom_list_query_enrichr_response,
)


class CustomListQueryEnrichrResultView(APIView):
    """
    Return the full Enrichr result for a successful
    directional CustomListQueryTask.

    Query params:
        taskUUID:
            CustomListQueryTask UUID.

        direction:
            up or down.

    File naming:
        {task_name}_mRNA_{direction}_enrichr.csv

    The endpoint returns all records. Ranking, searching and
    top-N selection are handled by the frontend.
    """

    def get(self, request, *args, **kwargs):
        try:
            task_uuid = get_required_task_uuid(
                request
            )

            task = get_task_or_raise(
                model_class=CustomListQueryTask,
                task_uuid=task_uuid,
                task_label="CustomListQueryTask",
            )

            require_success_task(
                task=task,
                task_label="Custom list query task",
            )

            if not bool(
                getattr(
                    task,
                    "has_mrna_direction",
                    False,
                )
            ):
                return Response(
                    {
                        "detail": (
                            "Directional Enrichr results are only "
                            "available when has_mrna_direction is true."
                        )
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            direction = normalize_enrichr_direction(
                request.query_params.get(
                    "direction"
                )
            )

            file_path, dataframe = (
                read_custom_list_query_enrichr_file(
                    task=task,
                    direction=direction,
                )
            )

            response_data = (
                build_custom_list_query_enrichr_response(
                    task=task,
                    direction=direction,
                    file_path=file_path,
                    dataframe=dataframe,
                )
            )

            return Response(
                response_data,
                status=status.HTTP_200_OK,
            )

        except WorkflowNetworkViewError as exc:
            return Response(
                {
                    "detail": exc.msg,
                },
                status=exc.status_code,
            )

        except (
            CustomListQueryEnrichrInputError,
            CustomListQueryEnrichrPathError,
        ) as exc:
            return Response(
                {
                    "detail": str(exc),
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        except FileNotFoundError as exc:
            return Response(
                {
                    "detail": str(exc),
                },
                status=status.HTTP_404_NOT_FOUND,
            )

        except Exception as exc:
            print(traceback.format_exc())

            return Response(
                {
                    "detail": (
                        f"Server error: {str(exc)}"
                    ),
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )
