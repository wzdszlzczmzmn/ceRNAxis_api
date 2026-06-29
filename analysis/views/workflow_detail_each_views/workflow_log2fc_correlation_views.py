import traceback

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.models import (
    PairedCohortTask,
    HybridReferenceTask,
)

from analysis.utils.workflow_detail_utils.workflow_network_view_utils import (
    WorkflowNetworkViewError,
    get_required_task_uuid,
    get_task_or_raise,
    require_success_task,
)

from analysis.utils.workflow_detail_utils.workflow_log2fc_background_utils import (
    PAIRED_COHORT_VALID_BACKGROUND_TYPES,
    HYBRID_REFERENCE_VALID_BACKGROUND_TYPES,
    WORKFLOW_LOG2FC_X_COL,
    WORKFLOW_LOG2FC_Y_COL,
    WorkflowLog2FCBackgroundPathError,
    WorkflowLog2FCBackgroundInputError,
    read_workflow_log2fc_background_file,
    get_workflow_available_background_types,
    build_workflow_log2fc_correlation_response_data,
)


class WorkflowLog2FCCorrelationBaseView(APIView):
    """
    Return core data for workflow log2FC background correlation plot.

    Query params:
        taskUUID: Workflow task UUID
        type: optional interaction type

    Input filename:
        {task_name}_ceRNA_background.csv

    Plot mapping:
        x = ceRNA_log2FC
        y = miRNA_log2FC
    """

    task_model = None
    task_type = None
    task_label = "Workflow task"

    valid_background_types = []
    x_col = WORKFLOW_LOG2FC_X_COL
    y_col = WORKFLOW_LOG2FC_Y_COL

    def get_extra_response_data(self, task) -> dict:
        return {}

    def get(self, request):
        try:
            if self.task_model is None:
                raise RuntimeError("Missing task_model.")

            if not self.task_type:
                raise RuntimeError("Missing task_type.")

            if not self.valid_background_types:
                raise RuntimeError("Missing valid_background_types.")

            task_uuid = get_required_task_uuid(request)

            requested_type = request.query_params.get("type", None)

            type_value = (
                str(requested_type).strip()
                if requested_type is not None
                else ""
            )

            task = get_task_or_raise(
                model_class=self.task_model,
                task_uuid=task_uuid,
                task_label=self.task_type,
            )

            require_success_task(
                task=task,
                task_label=self.task_label,
            )

            try:
                background_file, df = read_workflow_log2fc_background_file(task)
            except WorkflowLog2FCBackgroundPathError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_400_BAD_REQUEST,
                )
            except FileNotFoundError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_404_NOT_FOUND,
                )
            except WorkflowLog2FCBackgroundInputError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            available_types = get_workflow_available_background_types(
                df=df,
                valid_types=self.valid_background_types,
            )

            if not available_types:
                return Response(
                    {
                        "detail": (
                            "No supported background interaction type found "
                            "in ceRNA background file."
                        )
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not type_value:
                type_value = available_types[0]

            if type_value not in available_types:
                return Response(
                    {
                        "detail": (
                            "Invalid type for this task. Allowed values are: "
                            f"{', '.join(available_types)}."
                        ),
                        "available_types": available_types,
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            result = build_workflow_log2fc_correlation_response_data(
                task=task,
                task_type=self.task_type,
                df=df,
                type_value=type_value,
                background_file_name=background_file.name,
                available_types=available_types,
                x_col=self.x_col,
                y_col=self.y_col,
            )

            result.update(
                self.get_extra_response_data(task)
            )

            return Response(
                result,
                status=status.HTTP_200_OK,
            )

        except WorkflowNetworkViewError as e:
            return Response(
                {
                    "detail": e.msg,
                },
                status=e.status_code,
            )

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {
                    "detail": f"Server error: {str(e)}",
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )


class PairedCohortLog2FCCorrelationView(WorkflowLog2FCCorrelationBaseView):
    task_model = PairedCohortTask
    task_type = "PairedCohortTask"
    task_label = "Paired cohort task"

    valid_background_types = PAIRED_COHORT_VALID_BACKGROUND_TYPES

    def get_extra_response_data(self, task) -> dict:
        return {
            "map_info": task.map_info,
            "deg_method": task.deg_method,
            "use_padj": getattr(task, "use_padj", True),
        }


class HybridReferenceLog2FCCorrelationView(WorkflowLog2FCCorrelationBaseView):
    task_model = HybridReferenceTask
    task_type = "HybridReferenceTask"
    task_label = "Hybrid reference task"

    valid_background_types = HYBRID_REFERENCE_VALID_BACKGROUND_TYPES

    def get_extra_response_data(self, task) -> dict:
        return {
            "map_info": task.map_info,
            "tcga_type": task.tcga_type,
            "lncrna_type": task.lncrna_type,
            "deg_method": task.deg_method,
            "use_padj": getattr(task, "use_padj", True),
        }
