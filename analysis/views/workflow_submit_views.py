import traceback

from django.utils import timezone
from rest_framework.parsers import JSONParser, MultiPartParser, FormParser
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.models import CustomListQueryTask, PairedCohortTask, HybridReferenceTask
from analysis.slurm_sbtch import sbatch_custom_list_query_task, sbatch_paired_cohort_task, sbatch_hybrid_reference_task
from analysis.utils.custom_list_query_task_utils import normalize_rnas, CustomListQueryTaskInputError, \
    prepare_custom_list_query_workspace, CustomListQueryPathError
from analysis.utils.hybrid_reference_task_utils import validate_hybrid_reference_task_params, \
    HybridReferenceTaskInputError, HybridReferenceTaskPathError, HYBRID_REFERENCE_ALLOWED_FILE_FIELDS, \
    prepare_hybrid_reference_workspace, save_hybrid_reference_uploaded_input_files, \
    validate_hybrid_reference_file_contents
from analysis.utils.immune_annotation_path_utils import (
    ImmuneAnnotationPathError,
    validate_immune_annotation_file,
)
from analysis.utils.paired_cohort_task_utils import PAIRED_COHORT_ALLOWED_FILE_FIELDS, \
    prepare_paired_cohort_workspace, save_paired_cohort_uploaded_input_files, PairedCohortTaskInputError, \
    PairedCohortTaskPathError, validate_paired_cohort_file_contents, PAIRED_COHORT_REQUIRED_FILE_FIELDS


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
    PAIRED_COHORT_ALLOWED_CANCER_TYPES = [
        "MEL",
        "LUAD",
        "OS",
        "STAD",
        "BRCA",
        "CRC",
        "NSCLC",
        "OV",
        "LUSC",
        "UCEC",
        "CESC",
        "HCC",
        "AML",
        "ALL",
        "PRAD",
        "SCLC",
        "NBL",
        "MM",
        "Lymphoma",
        "PAAD",
        "",
    ]

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
                        "msg": (
                            "Invalid field: deg_method. "
                            "Allowed values: limma, deseq2."
                        ),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            cancer_type = str(request.data.get("cancer_type", "")).strip()
            use_padj = self.parse_bool_field(
                request,
                "use_padj",
                default=True,
            )

            if cancer_type not in self.PAIRED_COHORT_ALLOWED_CANCER_TYPES:
                return Response(
                    {
                        "success": False,
                        "msg": (
                            "Invalid field: cancer_type. "
                            "Allowed values are: "
                            f"{', '.join([x for x in self.PAIRED_COHORT_ALLOWED_CANCER_TYPES if x])} "
                            "or empty string."
                        ),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            missing_required_files = [
                field_name
                for field_name in PAIRED_COHORT_REQUIRED_FILE_FIELDS
                if field_name not in request.FILES
            ]

            if missing_required_files:
                return Response(
                    {
                        "success": False,
                        "msg": (
                            "Missing uploaded file(s): "
                            f"{', '.join(missing_required_files)}."
                        ),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if (
                    "lncrna_file" not in request.FILES
                    and "circrna_file" not in request.FILES
            ):
                return Response(
                    {
                        "success": False,
                        "msg": (
                            "At least one of lncrna_file or circrna_file is required."
                        ),
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
                logfc_cutoff_circrna = self.parse_float_field(
                    request,
                    "logfc_cutoff_circrna",
                )
                padj_cutoff_circrna = self.parse_float_field(
                    request,
                    "padj_cutoff_circrna",
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
                cancer_type=cancer_type,
                use_padj=use_padj,
                logfc_cutoff_mrna=logfc_cutoff_mrna,
                padj_cutoff_mrna=padj_cutoff_mrna,
                logfc_cutoff_mirna=logfc_cutoff_mirna,
                padj_cutoff_mirna=padj_cutoff_mirna,
                logfc_cutoff_lncrna=logfc_cutoff_lncrna,
                padj_cutoff_lncrna=padj_cutoff_lncrna,
                logfc_cutoff_circrna=logfc_cutoff_circrna,
                padj_cutoff_circrna=padj_cutoff_circrna,
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
                task.lncrna_file = saved_files.get("lncrna_file", "")
                task.circrna_file = saved_files.get("circrna_file", "")
                task.meta_file = saved_files["meta_file"]

                task.save(
                    update_fields=[
                        "mrna_file",
                        "mirna_file",
                        "lncrna_file",
                        "circrna_file",
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
                        "msg": (
                            f"Failed to prepare task workspace: {str(e)}"
                        ),
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
                        "cancer_type": task.cancer_type,
                        "use_padj": task.use_padj,
                        "files": {
                            "mrna_file": task.mrna_file,
                            "mirna_file": task.mirna_file,
                            "lncrna_file": task.lncrna_file,
                            "circrna_file": task.circrna_file,
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
                            "circRNA": {
                                "logfc_cutoff": task.logfc_cutoff_circrna,
                                "padj_cutoff": task.padj_cutoff_circrna,
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

    @staticmethod
    def parse_bool_field(request, field_name, default=None):
        raw_value = request.data.get(field_name, None)

        if raw_value is None or str(raw_value).strip() == "":
            if default is not None:
                return default
            raise ValueError(f"Missing field: {field_name}.")

        value = str(raw_value).strip().lower()

        if value in ["true", "1", "yes", "y"]:
            return True

        if value in ["false", "0", "no", "n"]:
            return False

        raise ValueError(
            f"Invalid boolean field: {field_name}. "
            "Allowed values: true, false."
        )


class HybridReferenceTaskSubmitView(APIView):
    parser_classes = [MultiPartParser, FormParser]

    def post(self, request):
        try:
            task_name = str(request.data.get("task_name", "")).strip()
            map_info = str(request.data.get("map_info", "")).strip()
            deg_method = str(request.data.get("deg_method", "")).strip()
            tcga_type = str(request.data.get("tcga_type", "")).strip()
            lncrna_type = str(request.data.get("lncrna_type", "")).strip()

            try:
                use_padj = self.parse_bool_field(request, "use_padj")
            except ValueError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

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

            if not tcga_type:
                return Response(
                    {
                        "success": False,
                        "msg": "Missing field: tcga_type.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not lncrna_type:
                return Response(
                    {
                        "success": False,
                        "msg": "Missing field: lncrna_type.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not deg_method:
                return Response(
                    {
                        "success": False,
                        "msg": "Missing field: deg_method.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                validate_hybrid_reference_task_params(
                    task_name=task_name,
                    tcga_type=tcga_type,
                    lncrna_type=lncrna_type,
                    deg_method=deg_method,
                    use_padj=use_padj,
                )
            except (
                    HybridReferenceTaskInputError,
                    HybridReferenceTaskPathError,
            ) as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            missing_files = [
                field_name
                for field_name in HYBRID_REFERENCE_ALLOWED_FILE_FIELDS
                if field_name not in request.FILES
            ]

            if missing_files:
                return Response(
                    {
                        "success": False,
                        "msg": (
                            "Missing uploaded file(s): "
                            f"{', '.join(missing_files)}."
                        ),
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

            except ValueError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            task = HybridReferenceTask.objects.create(
                user=request.data.get("user", ""),
                task_name=task_name,
                status=HybridReferenceTask.Status.Pending,
                map_info=map_info,
                tcga_type=tcga_type,
                lncrna_type=lncrna_type,
                use_padj=use_padj,
                deg_method=deg_method,
                logfc_cutoff_mrna=logfc_cutoff_mrna,
                padj_cutoff_mrna=padj_cutoff_mrna,
            )

            try:
                prepare_hybrid_reference_workspace(task)

                saved_files = save_hybrid_reference_uploaded_input_files(
                    task=task,
                    files=request.FILES,
                )

                validate_hybrid_reference_file_contents(task)

                task.mrna_file = saved_files["mrna_file"]
                task.meta_file = saved_files["meta_file"]

                task.save(
                    update_fields=[
                        "mrna_file",
                        "meta_file",
                    ]
                )

            except HybridReferenceTaskInputError as e:
                task.status = HybridReferenceTask.Status.Failed
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
                    HybridReferenceTaskPathError,
                    FileNotFoundError,
            ) as e:
                task.status = HybridReferenceTask.Status.Failed
                task.finish_time = timezone.now()
                task.save(update_fields=["status", "finish_time"])

                return Response(
                    {
                        "success": False,
                        "msg": (
                            f"Failed to prepare task workspace: {str(e)}"
                        ),
                    },
                    status=status.HTTP_500_INTERNAL_SERVER_ERROR,
                )

            submit_result = sbatch_hybrid_reference_task(task.uuid)

            if not submit_result["success"]:
                task.status = HybridReferenceTask.Status.Failed
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
                        "tcga_type": task.tcga_type,
                        "lncrna_type": task.lncrna_type,
                        "use_padj": task.use_padj,
                        "deg_method": task.deg_method,
                        "files": {
                            "mrna_file": task.mrna_file,
                            "meta_file": task.meta_file,
                        },
                        "cutoffs": {
                            "mRNA": {
                                "logfc_cutoff": task.logfc_cutoff_mrna,
                                "pvalue_cutoff": task.padj_cutoff_mrna,
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

    @staticmethod
    def parse_bool_field(request, field_name):
        raw_value = request.data.get(field_name, None)

        if raw_value is None or str(raw_value).strip() == "":
            raise ValueError(f"Missing field: {field_name}.")

        normalized_value = str(raw_value).strip().lower()

        if normalized_value in {"true", "1", "yes", "y"}:
            return True

        if normalized_value in {"false", "0", "no", "n"}:
            return False

        raise ValueError(
            f"Invalid boolean field: {field_name}. "
            "Allowed values are true or false."
        )
