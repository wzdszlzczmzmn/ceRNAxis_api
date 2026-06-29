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

from analysis.utils.workflow_detail_utils.workflow_axis_final_utils import (
    WorkflowAxisFinalPathError,
    WorkflowAxisFinalInputError,
    read_workflow_axis_final_file,
    normalize_workflow_axis_final_dataframe,
    serialize_workflow_axis_final_dataframe, HYBRID_REFERENCE_AXIS_FINAL_COLUMNS, PAIRED_COHORT_AXIS_FINAL_COLUMNS,
    WORKFLOW_AXIS_FINAL_NUMERIC_COLUMNS, build_axis_final_response_from_dataframe,
    PAIRED_COHORT_AXIS_FINAL_REQUIRED_COLUMNS, HYBRID_REFERENCE_AXIS_FINAL_REQUIRED_COLUMNS,
)


class WorkflowAxisFinalDataBaseView(APIView):
    """
    Return all rows from workflow ceRNA axis final result file.

    Query params:
        taskUUID: Workflow task UUID

    Input filename:
        {task_name}_ceRNA_axis_final.csv
    """

    task_model = None
    task_type = None
    task_label = "Workflow task"

    axis_final_columns = []
    axis_final_required_columns = []
    axis_final_numeric_columns = WORKFLOW_AXIS_FINAL_NUMERIC_COLUMNS

    def get_extra_response_data(self, task) -> dict:
        return {}

    def get(self, request):
        try:
            if self.task_model is None:
                raise RuntimeError("Missing task_model.")

            if not self.task_type:
                raise RuntimeError("Missing task_type.")

            if not self.axis_final_columns:
                raise RuntimeError("Missing axis_final_columns.")

            task_uuid = get_required_task_uuid(request)

            task = get_task_or_raise(
                model_class=self.task_model,
                task_uuid=task_uuid,
                task_label=self.task_type,
            )

            require_success_task(
                task=task,
                task_label=self.task_label,
            )

            required_columns = (
                self.axis_final_required_columns
                or self.axis_final_columns
            )

            try:
                axis_file, df = read_workflow_axis_final_file(
                    task=task,
                    required_columns=required_columns,
                )
            except WorkflowAxisFinalPathError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_400_BAD_REQUEST,
                )
            except FileNotFoundError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_404_NOT_FOUND,
                )
            except WorkflowAxisFinalInputError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            base_response = {
                "uuid": str(task.uuid),
                "task_type": self.task_type,
                "task_name": task.task_name,
            }

            base_response.update(
                self.get_extra_response_data(task)
            )

            response_data = build_axis_final_response_from_dataframe(
                df=df,
                axis_file_name=axis_file.name,
                columns=self.axis_final_columns,
                required_columns=required_columns,
                numeric_columns=self.axis_final_numeric_columns,
                base_response=base_response,
            )

            return Response(
                response_data,
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


class PairedCohortAxisFinalDataView(WorkflowAxisFinalDataBaseView):
    task_model = PairedCohortTask
    task_type = "PairedCohortTask"
    task_label = "Paired cohort task"

    axis_final_columns = PAIRED_COHORT_AXIS_FINAL_COLUMNS
    axis_final_required_columns = PAIRED_COHORT_AXIS_FINAL_REQUIRED_COLUMNS

    def get_extra_response_data(self, task) -> dict:
        return {
            "map_info": task.map_info,
            "deg_method": task.deg_method,
            "use_padj": getattr(task, "use_padj", True),
        }


class HybridReferenceAxisFinalDataView(WorkflowAxisFinalDataBaseView):
    task_model = HybridReferenceTask
    task_type = "HybridReferenceTask"
    task_label = "Hybrid reference task"

    axis_final_columns = HYBRID_REFERENCE_AXIS_FINAL_COLUMNS
    axis_final_required_columns = HYBRID_REFERENCE_AXIS_FINAL_REQUIRED_COLUMNS

    def get_extra_response_data(self, task) -> dict:
        return {
            "map_info": task.map_info,
            "tcga_type": task.tcga_type,
            "lncrna_type": task.lncrna_type,
            "deg_method": task.deg_method,
            "use_padj": getattr(task, "use_padj", True),
        }
