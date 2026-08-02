import traceback

from rest_framework import status
from rest_framework.response import Response
from rest_framework.views import APIView

from analysis.models import SCSTHybridReferenceTask
from analysis.utils.workflow_detail_utils.workflow_deg_volcano_utils import get_scst_workflow_deg_availability, \
    SCST_WORKFLOW_DEG_RNA_TYPES, SCST_WORKFLOW_DEG_SCOPES
from analysis.utils.workflow_detail_utils.workflow_log2fc_background_utils import \
    get_scst_workflow_background_availability, SCST_HYBRID_REFERENCE_VALID_BACKGROUND_TYPES

from analysis.utils.workflow_detail_utils.workflow_network_view_utils import (
    WorkflowNetworkViewError,
    get_required_task_uuid,
    get_task_or_raise,
    require_success_task,
)

from analysis.utils.workflow_detail_utils.workflow_viz_info_utils import (
    WorkflowVizInfoInputError,
    WorkflowVizInfoPathError,
    get_scst_hybrid_reference_group_info,
)


class SCSTHybridReferenceVizInfoView(APIView):
    """
    Return visualization options for an SC/ST Hybrid Reference task.

    Currently returned information:
    - group_col
    - unique group values
    - record count for each group
    - total group count
    - total sample/spot count

    Query parameters:
        taskUUID:
            SCSTHybridReferenceTask UUID

    Example:
        GET /api/analysis/scst_hybrid_reference_viz_info/
            ?taskUUID=df188272-9496-410a-8f28-ebef7bcdb672
    """

    task_model = SCSTHybridReferenceTask
    task_type = "SCSTHybridReferenceTask"
    task_label = "SC/ST hybrid reference task"

    def get(self, request):
        try:
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
                group_info = get_scst_hybrid_reference_group_info(
                    task
                )

                background_availability = (
                    get_scst_workflow_background_availability(
                        task=task,
                        group_values=group_info["group_values"],
                        valid_types=(
                            SCST_HYBRID_REFERENCE_VALID_BACKGROUND_TYPES
                        ),
                    )
                )

                deg_availability = (
                    get_scst_workflow_deg_availability(
                        task=task,
                        group_values=group_info["group_values"],
                        valid_rna_types=(
                            SCST_WORKFLOW_DEG_RNA_TYPES
                        ),
                        valid_scopes=(
                            SCST_WORKFLOW_DEG_SCOPES
                        ),
                    )
                )

                background_by_group = {
                    item.get("group_value"): item
                    for item in background_availability
                    if item.get("group_value")
                }

                deg_by_group = {
                    item.get("group_value"): item
                    for item in deg_availability
                    if item.get("group_value")
                }

                enriched_group_options = []

                for group_option in group_info["group_options"]:
                    group_value = group_option["value"]

                    background_info = background_by_group.get(
                        group_value,
                        {},
                    )

                    deg_info = deg_by_group.get(
                        group_value,
                        {},
                    )

                    enriched_group_options.append(
                        {
                            **group_option,

                            "background_file": background_info.get(
                                "background_file"
                            ),
                            "background_file_exists": (
                                background_info.get(
                                    "background_file_exists",
                                    False,
                                )
                            ),
                            "background_available": (
                                background_info.get(
                                    "background_available",
                                    False,
                                )
                            ),
                            "available_background_types": (
                                background_info.get(
                                    "available_background_types",
                                    [],
                                )
                            ),

                            "deg_available": deg_info.get(
                                "deg_available",
                                False,
                            ),
                            "available_deg_rna_types": (
                                deg_info.get(
                                    "available_deg_rna_types",
                                    [],
                                )
                            ),
                            "available_deg_scopes": (
                                deg_info.get(
                                    "available_deg_scopes",
                                    [],
                                )
                            ),
                            "deg_files": deg_info.get(
                                "deg_files",
                                {},
                            ),
                        }
                    )

                group_info["group_options"] = enriched_group_options

                group_info[
                    "available_deg_rna_types_by_group"
                ] = {
                    item.get("group_value"): item.get(
                        "available_deg_rna_types",
                        [],
                    )
                    for item in deg_availability
                    if item.get("group_value")
                }

                group_info[
                    "available_deg_scopes_by_group"
                ] = {
                    item.get("group_value"): item.get(
                        "available_deg_scopes",
                        [],
                    )
                    for item in deg_availability
                    if item.get("group_value")
                }

            except WorkflowVizInfoPathError as exc:
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

            except WorkflowVizInfoInputError as exc:
                return Response(
                    {
                        "detail": str(exc),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            response_data = {
                "uuid": str(task.uuid),
                "task_type": self.task_type,
                "task_name": task.task_name,

                "status": task.status,
                "status_label": task.get_status_display(),

                "data_type": task.data_type,
                "data_type_label": (
                    task.get_data_type_display()
                ),

                "tcga_type": task.tcga_type,
                "lncrna_type": task.lncrna_type,

                **group_info,
            }

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
