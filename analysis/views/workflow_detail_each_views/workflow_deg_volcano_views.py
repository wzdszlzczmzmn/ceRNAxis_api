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

from analysis.utils.workflow_detail_utils.workflow_deg_volcano_utils import (
    WORKFLOW_DEG_SCOPE_ALL,
    WORKFLOW_DEG_PAIRED_COHORT_RNA_TYPES,
    WORKFLOW_DEG_HYBRID_REFERENCE_RNA_TYPES,
    WORKFLOW_DEG_PAIRED_COHORT_SCOPES,
    WORKFLOW_DEG_HYBRID_REFERENCE_SCOPES,
    WorkflowDEGVolcanoPathError,
    WorkflowDEGVolcanoInputError,
    read_workflow_deg_file,
    build_workflow_deg_volcano_response_data,
    get_workflow_available_deg_rna_types,
    get_workflow_available_deg_scopes,
)


class WorkflowDEGVolcanoBaseView(APIView):
    """
    Return DEG volcano plot data for workflow task.

    Query params:
        taskUUID: Workflow task UUID
        rna_type: RNA type
        deg_scope: optional. default = all

    Default DEG filename:
        {task_name}_{deg_method}_{rna_type}.csv

    Intersect DEG filename:
        {task_name}_{deg_method}_{rna_type}_intersect.csv

    Required DEG columns:
        gene_name, log2FC, regulation, and one p-value column

    P-value source:
        task.use_padj == True  -> padj
        task.use_padj == False -> pvalue

    Response keeps unified field names:
        pvalue, neg_log10_pvalue
    """

    task_model = None
    task_type = None
    task_label = "Workflow task"

    valid_rna_types = []
    valid_deg_scopes = [
        WORKFLOW_DEG_SCOPE_ALL,
    ]

    default_rna_type = None
    default_deg_scope = WORKFLOW_DEG_SCOPE_ALL

    def get_extra_response_data(self, task) -> dict:
        return {}

    def get_default_rna_type(self) -> str | None:
        return self.default_rna_type

    def get_default_deg_scope(self) -> str:
        return self.default_deg_scope

    def get_use_padj(self, task) -> bool:
        return bool(getattr(task, "use_padj", True))

    def get(self, request):
        try:
            if self.task_model is None:
                raise RuntimeError("Missing task_model.")

            if not self.task_type:
                raise RuntimeError("Missing task_type.")

            if not self.valid_rna_types:
                raise RuntimeError("Missing valid_rna_types.")

            task_uuid = get_required_task_uuid(request)

            rna_type = str(
                request.query_params.get(
                    "rna_type",
                    self.get_default_rna_type() or "",
                )
            ).strip()

            deg_scope = str(
                request.query_params.get(
                    "deg_scope",
                    self.get_default_deg_scope(),
                )
            ).strip()

            if not rna_type:
                return Response(
                    {
                        "detail": "Missing query parameter: rna_type."
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if rna_type not in self.valid_rna_types:
                return Response(
                    {
                        "detail": (
                            "Invalid rna_type. Allowed values are: "
                            f"{', '.join(self.valid_rna_types)}."
                        ),
                        "valid_rna_types": self.valid_rna_types,
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if deg_scope not in self.valid_deg_scopes:
                return Response(
                    {
                        "detail": (
                            "Invalid deg_scope. Allowed values are: "
                            f"{', '.join(self.valid_deg_scopes)}."
                        ),
                        "valid_deg_scopes": self.valid_deg_scopes,
                    },
                    status=status.HTTP_400_BAD_REQUEST,
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

            available_deg_rna_types = get_workflow_available_deg_rna_types(
                task=task,
                valid_rna_types=self.valid_rna_types,
                deg_scope=WORKFLOW_DEG_SCOPE_ALL,
            )

            available_deg_scopes = get_workflow_available_deg_scopes(
                task=task,
                rna_type=rna_type,
                valid_scopes=self.valid_deg_scopes,
            )

            try:
                deg_file, df = read_workflow_deg_file(
                    task=task,
                    rna_type=rna_type,
                    deg_scope=deg_scope,
                )
            except WorkflowDEGVolcanoPathError as e:
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
                        "available_deg_rna_types": available_deg_rna_types,
                        "available_deg_scopes": available_deg_scopes,
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )
            except WorkflowDEGVolcanoInputError as e:
                return Response(
                    {
                        "detail": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            use_padj = self.get_use_padj(task)

            try:
                response_data = build_workflow_deg_volcano_response_data(
                    task=task,
                    task_type=self.task_type,
                    rna_type=rna_type,
                    deg_scope=deg_scope,
                    deg_file_name=deg_file.name,
                    df=df,
                    use_padj=use_padj,
                )
            except WorkflowDEGVolcanoInputError as e:
                return Response(
                    {
                        "detail": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            response_data.update(
                {
                    "available_deg_rna_types": available_deg_rna_types,
                    "available_deg_scopes": available_deg_scopes,
                }
            )

            response_data.update(
                self.get_extra_response_data(task)
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


class PairedCohortDEGVolcanoView(WorkflowDEGVolcanoBaseView):
    task_model = PairedCohortTask
    task_type = "PairedCohortTask"
    task_label = "Paired cohort task"

    valid_rna_types = WORKFLOW_DEG_PAIRED_COHORT_RNA_TYPES
    valid_deg_scopes = WORKFLOW_DEG_PAIRED_COHORT_SCOPES

    default_rna_type = "mRNA"
    default_deg_scope = WORKFLOW_DEG_SCOPE_ALL

    def get_extra_response_data(self, task) -> dict:
        return {
            "map_info": task.map_info,
            "deg_method": task.deg_method,
            "use_padj": getattr(task, "use_padj", True),
        }


class HybridReferenceDEGVolcanoView(WorkflowDEGVolcanoBaseView):
    task_model = HybridReferenceTask
    task_type = "HybridReferenceTask"
    task_label = "Hybrid reference task"

    valid_rna_types = WORKFLOW_DEG_HYBRID_REFERENCE_RNA_TYPES
    valid_deg_scopes = WORKFLOW_DEG_HYBRID_REFERENCE_SCOPES

    default_rna_type = "mRNA"
    default_deg_scope = WORKFLOW_DEG_SCOPE_ALL

    def get_extra_response_data(self, task) -> dict:
        return {
            "map_info": task.map_info,
            "tcga_type": task.tcga_type,
            "lncrna_type": task.lncrna_type,
            "deg_method": task.deg_method,
            "use_padj": getattr(task, "use_padj", True),
        }
