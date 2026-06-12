import traceback

from django.utils import timezone
from rest_framework.parsers import JSONParser
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.models import CustomListQueryTask
from analysis.slurm_sbtch import sbatch_custom_list_query_task
from analysis.utils.custom_list_query_task_utils import normalize_rnas, CustomListQueryTaskInputError, \
    prepare_custom_list_query_workspace, CustomListQueryPathError
from analysis.utils.immune_annotation_path_utils import (
    ImmuneAnnotationPathError,
    validate_immune_annotation_file,
)


class CustomListQueryTaskSubmitView(APIView):
    parser_classes = [JSONParser]

    def post(self, request):
        try:
            task_name = str(request.data.get("task_name", "")).strip()
            map_info = str(request.data.get("map_info", "")).strip()
            rnas = request.data.get("rnas")

            if not task_name:
                return Response(
                    {
                        "success": False,
                        "msg": "Missing field: task_name.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not map_info:
                return Response(
                    {
                        "success": False,
                        "msg": "Missing field: map_info.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                validate_immune_annotation_file(map_info)
            except ImmuneAnnotationPathError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )
            except FileNotFoundError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            try:
                normalized_rnas = normalize_rnas(rnas)
            except CustomListQueryTaskInputError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            task = CustomListQueryTask.objects.create(
                user=request.data.get("user", ""),
                task_name=task_name,
                status=CustomListQueryTask.Status.Pending,
                map_info=map_info,
                rnas=normalized_rnas,
            )

            try:
                workspace = prepare_custom_list_query_workspace(task)
            except (OSError, CustomListQueryTaskInputError, CustomListQueryPathError) as e:
                task.status = CustomListQueryTask.Status.Failed
                task.finish_time = timezone.now()
                task.save(update_fields=["status", "finish_time"])

                return Response(
                    {
                        "success": False,
                        "msg": f"Failed to prepare task workspace: {str(e)}",
                    },
                    status=status.HTTP_500_INTERNAL_SERVER_ERROR,
                )

            submit_result = sbatch_custom_list_query_task(task.uuid)

            if not submit_result["success"]:
                task.status = CustomListQueryTask.Status.Failed
                task.finish_time = timezone.now()
                task.save(update_fields=["status", "finish_time"])

                return Response(
                    {
                        "success": False,
                        "msg": submit_result["msg"],
                    },
                    status=status.HTTP_500_INTERNAL_SERVER_ERROR,
                )

            return Response(
                {
                    "success": True,
                    "msg": "Task submitted successfully.",
                    "data": {
                        "uuid": str(task.uuid),
                        "task_name": task.task_name,
                        "user": task.user,
                        "status": task.get_status_display(),
                        "create_time": timezone.localtime(
                            task.create_time
                        ).strftime("%Y-%m-%d %H:%M:%S"),
                        "map_info": task.map_info,
                        "rna_counts": {
                            "miRNA": task.miRNA_count,
                            "mRNA": task.mRNA_count,
                            "lncRNA": task.lncRNA_count,
                            "circRNA": task.circRNA_count,
                            "total": task.total_rna_count,
                        },
                    },
                },
                status=status.HTTP_201_CREATED,
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
