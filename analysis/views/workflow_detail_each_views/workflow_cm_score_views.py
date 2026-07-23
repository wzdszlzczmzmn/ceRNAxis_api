import traceback

from rest_framework import status
from rest_framework.response import Response
from rest_framework.views import APIView

from analysis.models import (
    CustomListQueryTask,
    PairedCohortTask,
    HybridReferenceTask,
)

from analysis.utils.workflow_detail_utils.workflow_network_view_utils import (
    WorkflowNetworkViewError,
    get_required_task_uuid,
    get_task_or_raise,
    require_success_task,
)

from analysis.utils.workflow_detail_utils.workflow_cm_score_utils import (
    WorkflowCMScoreInputError,
    WorkflowCMScorePathError,
    build_workflow_cm_score_options_response, validate_cm_score_item_value, read_workflow_cm_score_file,
    build_workflow_cm_score_result_response,
)


WORKFLOW_CM_SCORE_TASK_CONFIG = {
    "CustomListQueryTask": {
        "model": CustomListQueryTask,
        "label": "Custom list query task",
    },
    "PairedCohortTask": {
        "model": PairedCohortTask,
        "label": "Paired cohort task",
    },
    "HybridReferenceTask": {
        "model": HybridReferenceTask,
        "label": "Hybrid reference task",
    },
}


def get_required_task_type(request) -> str:
    task_type = str(
        request.query_params.get(
            "taskType",
            "",
        )
    ).strip()

    if not task_type:
        raise WorkflowCMScoreInputError(
            "Missing required query parameter: taskType."
        )

    if task_type not in WORKFLOW_CM_SCORE_TASK_CONFIG:
        raise WorkflowCMScoreInputError(
            "Unsupported taskType. Allowed values are: "
            f"{', '.join(WORKFLOW_CM_SCORE_TASK_CONFIG)}."
        )

    return task_type


class WorkflowCMScoreOptionsView(APIView):
    """
    Return available CM-score result items for a workflow task.

    This endpoint is shared by:
        - CustomListQueryTask
        - PairedCohortTask
        - HybridReferenceTask

    Query parameters:
        taskType:
            CustomListQueryTask
            PairedCohortTask
            HybridReferenceTask

        taskUUID:
            Workflow task UUID

    File layout:
        task_output_dir/
            CM_results/
                {item}_CM_scores.csv

    The value represented by `item` depends on the module:
        Module 1:
            gene

        Module 2 / Module 3:
            axis ID

    The API deliberately uses the generic term `item`.
    """

    def get(self, request, *args, **kwargs):
        try:
            task_type = get_required_task_type(
                request
            )

            task_uuid = get_required_task_uuid(
                request
            )

            task_config = (
                WORKFLOW_CM_SCORE_TASK_CONFIG[
                    task_type
                ]
            )

            task = get_task_or_raise(
                model_class=task_config["model"],
                task_uuid=task_uuid,
                task_label=task_type,
            )

            require_success_task(
                task=task,
                task_label=task_config["label"],
            )

            response_data = (
                build_workflow_cm_score_options_response(
                    task=task,
                    task_type=task_type,
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
            WorkflowCMScoreInputError,
            WorkflowCMScorePathError,
        ) as exc:
            return Response(
                {
                    "detail": str(exc),
                },
                status=status.HTTP_400_BAD_REQUEST,
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


class WorkflowCMScoreResultView(APIView):
    """
    Return all drug-signature rows from a selected CM-score file.

    Shared by:
        CustomListQueryTask
        PairedCohortTask
        HybridReferenceTask

    Query parameters:
        taskType:
            CustomListQueryTask
            PairedCohortTask
            HybridReferenceTask

        taskUUID:
            Workflow task UUID.

        item:
            Generic CM-score item identifier.

            Module 1:
                gene

            Module 2 / Module 3:
                axis ID

    File layout:
        task_output_dir/
            CM_results/
                {item}_CM_scores.csv
    """

    def get(self, request, *args, **kwargs):
        try:
            task_type = get_required_task_type(
                request
            )

            task_uuid = get_required_task_uuid(
                request
            )

            item_value = validate_cm_score_item_value(
                request.query_params.get("item")
            )

            task_config = (
                WORKFLOW_CM_SCORE_TASK_CONFIG[
                    task_type
                ]
            )

            task = get_task_or_raise(
                model_class=task_config["model"],
                task_uuid=task_uuid,
                task_label=task_type,
            )

            require_success_task(
                task=task,
                task_label=task_config["label"],
            )

            file_path, dataframe = (
                read_workflow_cm_score_file(
                    task=task,
                    item_value=item_value,
                )
            )

            response_data = (
                build_workflow_cm_score_result_response(
                    task=task,
                    task_type=task_type,
                    item_value=item_value,
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
            WorkflowCMScoreInputError,
            WorkflowCMScorePathError,
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
                    "detail": f"Server error: {str(exc)}",
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )
