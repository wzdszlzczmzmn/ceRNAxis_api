import traceback
import uuid as uuid_lib

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.models import PairedCohortTask, HybridReferenceTask, SCSTHybridReferenceTask
from analysis.utils.paired_cohort_task_utils import (
    PairedCohortTaskInputError,
    PairedCohortTaskPathError,
    build_paired_cohort_survival_km_data,
)
from analysis.utils.hybrid_reference_task_utils import (
    HybridReferenceTaskInputError,
    HybridReferenceTaskPathError,
    build_hybrid_reference_survival_km_data, SCSTHybridReferenceTaskPathError, SCSTHybridReferenceTaskInputError,
    build_scst_hybrid_reference_survival_km_data,
)
from analysis.utils.workflow_detail_utils.workflow_viz_info_utils import WorkflowVizInfoPathError, \
    WorkflowVizInfoInputError, validate_scst_hybrid_reference_h5ad_group_value


class BaseWorkflowSurvivalKMDataView(APIView):
    task_model = None
    task_label = "Workflow task"
    not_success_message = (
        "Task is not completed successfully."
    )

    build_result_func = None

    path_error_classes = ()
    input_error_classes = ()

    def get_task(
        self,
        task_uuid: str,
    ):
        return self.task_model.objects.get(
            uuid=task_uuid
        )

    def get_status_success_value(self):
        return self.task_model.Status.Success

    def build_result(
        self,
        *,
        request,
        task,
    ) -> dict:
        if self.build_result_func is None:
            raise RuntimeError(
                "Missing build_result_func."
            )

        return self.build_result_func(
            task
        )

    def get(self, request):
        try:
            task_uuid = str(
                request.query_params.get(
                    "taskUUID",
                    "",
                )
            ).strip()

            if not task_uuid:
                return Response(
                    {
                        "detail":
                            "Missing query parameter: taskUUID."
                    },
                    status=(
                        status.HTTP_400_BAD_REQUEST
                    ),
                )

            try:
                uuid_lib.UUID(
                    task_uuid
                )
            except ValueError:
                return Response(
                    {
                        "detail":
                            "Invalid Task UUID format."
                    },
                    status=(
                        status.HTTP_400_BAD_REQUEST
                    ),
                )

            try:
                task = self.get_task(
                    task_uuid
                )
            except self.task_model.DoesNotExist:
                return Response(
                    {
                        "detail": (
                            f"{self.task_label} "
                            f"not found: {task_uuid}."
                        )
                    },
                    status=(
                        status.HTTP_404_NOT_FOUND
                    ),
                )

            if (
                task.status
                != self.get_status_success_value()
            ):
                return Response(
                    {
                        "detail": (
                            f"{self.not_success_message} "
                            "Current status: "
                            f"{task.get_status_display()}."
                        )
                    },
                    status=(
                        status.HTTP_400_BAD_REQUEST
                    ),
                )

            try:
                result = self.build_result(
                    request=request,
                    task=task,
                )

            except self.path_error_classes as exc:
                return Response(
                    {
                        "detail": str(exc)
                    },
                    status=(
                        status.HTTP_400_BAD_REQUEST
                    ),
                )

            except FileNotFoundError as exc:
                return Response(
                    {
                        "detail": str(exc)
                    },
                    status=(
                        status.HTTP_404_NOT_FOUND
                    ),
                )

            except self.input_error_classes as exc:
                return Response(
                    {
                        "detail": str(exc)
                    },
                    status=(
                        status.HTTP_400_BAD_REQUEST
                    ),
                )

            return Response(
                result,
                status=status.HTTP_200_OK,
            )

        except Exception as exc:
            print(
                traceback.format_exc()
            )

            return Response(
                {
                    "detail":
                        f"Server error: {str(exc)}"
                },
                status=(
                    status.HTTP_500_INTERNAL_SERVER_ERROR
                ),
            )


class PairedCohortSurvivalKMDataView(BaseWorkflowSurvivalKMDataView):
    """
    Return Kaplan-Meier survival curve data for PairedCohortTask.

    Query params:
        taskUUID: PairedCohortTask UUID

    Input filename:
        {task_name}_survival_analysis.csv
    """

    task_model = PairedCohortTask
    task_label = "PairedCohortTask"
    not_success_message = "Paired cohort task is not completed successfully."
    build_result_func = staticmethod(build_paired_cohort_survival_km_data)
    path_error_classes = (PairedCohortTaskPathError,)
    input_error_classes = (PairedCohortTaskInputError,)


class HybridReferenceSurvivalKMDataView(BaseWorkflowSurvivalKMDataView):
    """
    Return Kaplan-Meier survival curve data for HybridReferenceTask.

    Query params:
        taskUUID: HybridReferenceTask UUID

    Input filename:
        {task_name}_survival_analysis.csv
    """

    task_model = HybridReferenceTask
    task_label = "HybridReferenceTask"
    not_success_message = "Hybrid reference task is not completed successfully."
    build_result_func = staticmethod(build_hybrid_reference_survival_km_data)
    path_error_classes = (HybridReferenceTaskPathError,)
    input_error_classes = (HybridReferenceTaskInputError,)


class SCSTHybridReferenceSurvivalKMDataView(
    BaseWorkflowSurvivalKMDataView
):
    """
    Return Kaplan-Meier survival curve data
    for SCSTHybridReferenceTask.

    Query params:
        taskUUID:
            SCSTHybridReferenceTask UUID

        groupValue:
            Selected SC/ST group value

    Input filename:
        {task_name}_survival_analysis_{groupValue}.csv
    """

    task_model = SCSTHybridReferenceTask
    task_label = "SCSTHybridReferenceTask"
    not_success_message = "SC/ST hybrid reference task is not completed successfully."
    path_error_classes = (
        SCSTHybridReferenceTaskPathError,
        WorkflowVizInfoPathError,
    )
    input_error_classes = (
        SCSTHybridReferenceTaskInputError,
        WorkflowVizInfoInputError,
    )

    def build_result(
        self,
        *,
        request,
        task,
    ) -> dict:
        group_value = str(
            request.query_params.get(
                "groupValue",
                "",
            )
        ).strip()

        if not group_value:
            raise WorkflowVizInfoInputError(
                "Missing query parameter: groupValue."
            )

        group_value = (
            validate_scst_hybrid_reference_h5ad_group_value(
                task=task,
                group_value=group_value,
            )
        )

        return (
            build_scst_hybrid_reference_survival_km_data(
                task=task,
                group_value=group_value,
            )
        )
