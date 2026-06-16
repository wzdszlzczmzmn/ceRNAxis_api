import traceback

from django.utils import timezone
from rest_framework.parsers import JSONParser, MultiPartParser, FormParser
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.models import CustomListQueryTask, PairedCohortTask
from analysis.slurm_sbtch import sbatch_custom_list_query_task, sbatch_paired_cohort_task
from analysis.utils.custom_list_query_task_utils import normalize_rnas, CustomListQueryTaskInputError, \
    prepare_custom_list_query_workspace, CustomListQueryPathError
from analysis.utils.immune_annotation_path_utils import (
    ImmuneAnnotationPathError,
    validate_immune_annotation_file,
)
from analysis.utils.paired_cohort_task_utils import PAIRED_COHORT_ALLOWED_FILE_FIELDS, \
    prepare_paired_cohort_workspace, save_paired_cohort_uploaded_input_files, PairedCohortTaskInputError, \
    PairedCohortTaskPathError, validate_paired_cohort_file_contents


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


class PairedCohortTaskSubmitView(APIView):
    parser_classes = [MultiPartParser, FormParser]

    def post(self, request):
        try:
            task_name = str(request.data.get("task_name", "")).strip()
            map_info = str(request.data.get("map_info", "")).strip()
            deg_method = str(request.data.get("deg_method", "")).strip()

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

            if not deg_method:
                return Response(
                    {
                        "success": False,
                        "msg": "Missing field: deg_method.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if deg_method not in ["limma", "deseq2"]:
                return Response(
                    {
                        "success": False,
                        "msg": "Invalid field: deg_method. Allowed values: limma, deseq2.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            missing_files = [
                field_name
                for field_name in PAIRED_COHORT_ALLOWED_FILE_FIELDS
                if field_name not in request.FILES
            ]

            if missing_files:
                return Response(
                    {
                        "success": False,
                        "msg": f"Missing uploaded file(s): {', '.join(missing_files)}.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                logfc_cutoff_mrna = self.parse_float_field(
                    request,
                    "logfc_cutoff_mrna",
                )
                padj_cutoff_mrna = self.parse_float_field(
                    request,
                    "padj_cutoff_mrna",
                )
                logfc_cutoff_mirna = self.parse_float_field(
                    request,
                    "logfc_cutoff_mirna",
                )
                padj_cutoff_mirna = self.parse_float_field(
                    request,
                    "padj_cutoff_mirna",
                )
                logfc_cutoff_lncrna = self.parse_float_field(
                    request,
                    "logfc_cutoff_lncrna",
                )
                padj_cutoff_lncrna = self.parse_float_field(
                    request,
                    "padj_cutoff_lncrna",
                )
            except ValueError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            task = PairedCohortTask.objects.create(
                user=request.data.get("user", ""),
                task_name=task_name,
                status=PairedCohortTask.Status.Pending,
                map_info=map_info,
                deg_method=deg_method,
                logfc_cutoff_mrna=logfc_cutoff_mrna,
                padj_cutoff_mrna=padj_cutoff_mrna,
                logfc_cutoff_mirna=logfc_cutoff_mirna,
                padj_cutoff_mirna=padj_cutoff_mirna,
                logfc_cutoff_lncrna=logfc_cutoff_lncrna,
                padj_cutoff_lncrna=padj_cutoff_lncrna,
            )

            try:
                prepare_paired_cohort_workspace(task)

                saved_files = save_paired_cohort_uploaded_input_files(
                    task=task,
                    files=request.FILES,
                )

                validate_paired_cohort_file_contents(task)

                task.mrna_file = saved_files["mrna_file"]
                task.mirna_file = saved_files["mirna_file"]
                task.lncrna_file = saved_files["lncrna_file"]
                task.meta_file = saved_files["meta_file"]
                task.save(
                    update_fields=[
                        "mrna_file",
                        "mirna_file",
                        "lncrna_file",
                        "meta_file",
                    ]
                )

            except PairedCohortTaskInputError as e:
                task.status = PairedCohortTask.Status.Failed
                task.finish_time = timezone.now()
                task.save(update_fields=["status", "finish_time"])

                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            except (
                    OSError,
                    PairedCohortTaskPathError,
                    FileNotFoundError,
            ) as e:
                task.status = PairedCohortTask.Status.Failed
                task.finish_time = timezone.now()
                task.save(update_fields=["status", "finish_time"])

                return Response(
                    {

                        "success": False,

                        "msg": f"Failed to prepare task workspace: {str(e)}",

                    },
                    status=status.HTTP_500_INTERNAL_SERVER_ERROR,
                )

            submit_result = sbatch_paired_cohort_task(task.uuid)

            if not submit_result["success"]:
                task.status = PairedCohortTask.Status.Failed
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
                        "deg_method": task.deg_method,
                        "files": {
                            "mrna_file": task.mrna_file,
                            "mirna_file": task.mirna_file,
                            "lncrna_file": task.lncrna_file,
                            "meta_file": task.meta_file,
                        },
                        "cutoffs": {
                            "mRNA": {
                                "logfc_cutoff": task.logfc_cutoff_mrna,
                                "padj_cutoff": task.padj_cutoff_mrna,
                            },
                            "miRNA": {
                                "logfc_cutoff": task.logfc_cutoff_mirna,
                                "padj_cutoff": task.padj_cutoff_mirna,
                            },
                            "lncRNA": {
                                "logfc_cutoff": task.logfc_cutoff_lncrna,
                                "padj_cutoff": task.padj_cutoff_lncrna,
                            },
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

    @staticmethod
    def parse_float_field(request, field_name):
        raw_value = request.data.get(field_name, None)

        if raw_value is None or str(raw_value).strip() == "":
            raise ValueError(f"Missing field: {field_name}.")

        try:
            value = float(raw_value)
        except (TypeError, ValueError):
            raise ValueError(f"Invalid numeric field: {field_name}.")

        if field_name.startswith("logfc_cutoff"):
            if value < 0:
                raise ValueError(
                    f"{field_name} must be greater than or equal to 0."
                )

        if field_name.startswith("padj_cutoff"):
            if value <= 0 or value > 1:
                raise ValueError(
                    f"{field_name} must be in the range (0, 1]."
                )

        return value
