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

from analysis.utils.workflow_detail_utils.workflow_cmap_utils import (
    WORKFLOW_CMAP_REQUIRED_COLUMNS,
    HYBRID_REFERENCE_CMAP_REQUIRED_COLUMNS,
    WorkflowCMapPathError,
    WorkflowCMapInputError,
    read_workflow_cmap_file,
    normalize_workflow_cmap_dataframe,
    serialize_workflow_cmap_dataframe, build_cmap_response_from_dataframe,
)


class WorkflowCMapResultBaseView(APIView):
    """
    Return full CMap result table for workflow task.

    Query params:
        taskUUID: Workflow task UUID

    CMap filename rule:
        {task_name}_CMap.csv
    """

    task_model = None
    task_type = None
    task_label = "Workflow task"
    required_columns = WORKFLOW_CMAP_REQUIRED_COLUMNS

    def get_extra_response_data(self, task) -> dict:
        return {}

    def get(self, request):
        try:
            if self.task_model is None:
                raise RuntimeError("Missing task_model.")

            if not self.task_type:
                raise RuntimeError("Missing task_type.")

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

            try:
                cmap_file, df = read_workflow_cmap_file(
                    task=task,
                    required_columns=self.required_columns,
                )
            except WorkflowCMapPathError as e:
                return Response(
                    {
                        "detail": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )
            except FileNotFoundError as e:
                return Response(
                    {
                        "detail": str(e),
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )
            except WorkflowCMapInputError as e:
                return Response(
                    {
                        "detail": str(e),
                    },
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

            response_data = build_cmap_response_from_dataframe(
                df=df,
                cmap_file_name=cmap_file.name,
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


class PairedCohortCMapResultView(WorkflowCMapResultBaseView):
    task_model = PairedCohortTask
    task_type = "PairedCohortTask"
    task_label = "Paired cohort task"

    # Paired Cohort CMap 当前不强制列结构，保持旧逻辑。
    required_columns = WORKFLOW_CMAP_REQUIRED_COLUMNS

    def get_extra_response_data(self, task) -> dict:
        return {
            "map_info": task.map_info,
            "deg_method": task.deg_method,
            "use_padj": getattr(task, "use_padj", True),
        }


class HybridReferenceCMapResultView(WorkflowCMapResultBaseView):
    task_model = HybridReferenceTask
    task_type = "HybridReferenceTask"
    task_label = "Hybrid reference task"

    required_columns = HYBRID_REFERENCE_CMAP_REQUIRED_COLUMNS

    def get_extra_response_data(self, task) -> dict:
        return {
            "map_info": task.map_info,
            "tcga_type": task.tcga_type,
            "lncrna_type": task.lncrna_type,
            "deg_method": task.deg_method,
            "use_padj": getattr(task, "use_padj", True),
        }
