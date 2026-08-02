import traceback

from rest_framework import status
from rest_framework.response import Response
from rest_framework.views import APIView

from analysis.models import (
    HybridReferenceTask,
    PairedCohortTask,
)

from analysis.utils.workflow_detail_utils.workflow_network_view_utils import (
    WorkflowNetworkViewError,
    get_required_task_uuid,
    get_task_or_raise,
    require_success_task,
)

from analysis.utils.workflow_detail_utils.workflow_sponge_utils import (
    WORKFLOW_SPONGE_COLUMNS,
    WORKFLOW_SPONGE_NUMERIC_COLUMNS,
    WORKFLOW_SPONGE_REQUIRED_COLUMNS,
    WorkflowSpongeInputError,
    WorkflowSpongePathError,
    build_sponge_response_from_dataframe,
    read_workflow_sponge_file,
)


class WorkflowSpongeDataBaseView(APIView):
    """
    Return all rows from a workflow Sponge result file.

    Query parameters:
        taskUUID: workflow task UUID

    Expected input filename:
        {task_name}_sponge_result.csv
    """

    task_model = None
    task_type = None
    task_label = "Workflow task"

    sponge_columns = WORKFLOW_SPONGE_COLUMNS

    sponge_required_columns = (
        WORKFLOW_SPONGE_REQUIRED_COLUMNS
    )

    sponge_numeric_columns = (
        WORKFLOW_SPONGE_NUMERIC_COLUMNS
    )

    def get_extra_response_data(self, task) -> dict:
        """
        Return task-type-specific response metadata.
        """
        return {}

    def should_include_project_matches(
        self,
        request,
    ) -> bool:
        """
        Include recurrent-reference matches by default.

        Use:
            ?include_project_matches=false
        to disable matching.
        """
        value = request.query_params.get(
            "include_project_matches"
        )

        if value is None:
            return True

        return str(value).strip().lower() not in {
            "0",
            "false",
            "no",
            "off",
        }

    def get_max_matches_per_axis(
        self,
        request,
    ) -> int | None:
        value = request.query_params.get(
            "max_matches_per_axis"
        )

        if value in {None, ""}:
            return None

        try:
            value = int(value)
        except (TypeError, ValueError):
            return None

        return value if value > 0 else None

    def get_project_match_scope(
        self,
        request,
    ) -> str:
        """
        Sponge results match Sponge reference contexts by default.

        Supported by the recurrent Axis service:
            sponge, both, any, axis_final
        """
        return str(
            request.query_params.get(
                "project_match_scope",
                "sponge",
            )
        ).strip().lower()

    def get(self, request):
        try:
            self.validate_view_configuration()

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
                self.sponge_required_columns
                or self.sponge_columns
            )

            try:
                sponge_file, df = (
                    read_workflow_sponge_file(
                        task=task,
                        required_columns=required_columns,
                    )
                )
            except WorkflowSpongePathError as exc:
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
            except WorkflowSpongeInputError as exc:
                return Response(
                    {
                        "detail": str(exc),
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

            response_data = (
                build_sponge_response_from_dataframe(
                    df=df,
                    sponge_file_name=sponge_file.name,
                    columns=self.sponge_columns,
                    required_columns=required_columns,
                    numeric_columns=(
                        self.sponge_numeric_columns
                    ),
                    base_response=base_response,
                    include_project_matches=(
                        self.should_include_project_matches(
                            request
                        )
                    ),
                    max_matches_per_axis=(
                        self.get_max_matches_per_axis(
                            request
                        )
                    ),
                    project_match_scope=(
                        self.get_project_match_scope(
                            request
                        )
                    ),
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

        except Exception as exc:
            print(traceback.format_exc())

            return Response(
                {
                    "detail": (
                        f"Server error: {str(exc)}"
                    ),
                },
                status=(
                    status.HTTP_500_INTERNAL_SERVER_ERROR
                ),
            )

    def validate_view_configuration(self) -> None:
        """
        Validate required subclass configuration.
        """
        if self.task_model is None:
            raise RuntimeError(
                "Missing task_model."
            )

        if not self.task_type:
            raise RuntimeError(
                "Missing task_type."
            )

        if not self.sponge_columns:
            raise RuntimeError(
                "Missing sponge_columns."
            )


class PairedCohortSpongeDataView(
    WorkflowSpongeDataBaseView
):
    task_model = PairedCohortTask
    task_type = "PairedCohortTask"
    task_label = "Paired cohort task"

    def get_extra_response_data(self, task) -> dict:
        return {
            "map_info": task.map_info,
            "deg_method": task.deg_method,
            "cancer_type": task.cancer_type,
            "use_padj": getattr(
                task,
                "use_padj",
                True,
            ),
        }


class HybridReferenceSpongeDataView(
    WorkflowSpongeDataBaseView
):
    task_model = HybridReferenceTask
    task_type = "HybridReferenceTask"
    task_label = "Hybrid reference task"

    def get_extra_response_data(self, task) -> dict:
        return {
            "map_info": task.map_info,
            "tcga_type": task.tcga_type,
            "lncrna_type": task.lncrna_type,
            "deg_method": task.deg_method,
            "use_padj": getattr(
                task,
                "use_padj",
                True,
            ),
        }
