import traceback

import pandas as pd
from django.http import FileResponse
from rest_framework.parsers import JSONParser
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from django.utils import timezone

from analysis.models import CustomListQueryTask, PairedCohortTask
from analysis.slurm_sbtch import sbatch_custom_list_query_task, sbatch_paired_cohort_task
from analysis.utils.custom_list_query_demo_utils import load_custom_list_query_demo_input, CustomListQueryDemoPathError, \
    CustomListQueryDemoConfigError
from analysis.utils.custom_list_query_task_utils import normalize_rnas, CustomListQueryTaskInputError, \
    prepare_custom_list_query_workspace, CustomListQueryPathError, validate_task_name_for_filename
from analysis.utils.immune_annotation_path_utils import validate_immune_annotation_file, ImmuneAnnotationPathError
from analysis.utils.paired_cohort_demo_utils import (
    PAIRED_COHORT_DEMO_VALID_RNA_TYPES,
    get_paired_cohort_demo_dir,
    load_paired_cohort_demo_manifest,
    get_paired_cohort_demo_meta_file_path,
    get_paired_cohort_demo_expression_file_path,
    validate_demo_file_exists,
    read_parquet_columns,
    PairedCohortDemoError,
    PairedCohortDemoPathError,
    PairedCohortDemoManifestError, validate_paired_cohort_demo_data_archive, PAIRED_COHORT_DEMO_DATA_ARCHIVE_NAME,
    load_paired_cohort_demo_input, PairedCohortDemoConfigError, get_paired_cohort_demo_cutoff_fields,
    copy_paired_cohort_demo_input_files_to_task,
)
from analysis.utils.paired_cohort_task_utils import PairedCohortTaskPathError, prepare_paired_cohort_workspace, \
    validate_paired_cohort_file_contents, PairedCohortTaskInputError

PAIRED_COHORT_DEMO_MAX_SELECTED_GENES = 30
DEFAULT_DEMO_EXPRESSION_FILE_FORMAT = "parquet"


class CustomListQueryDemoRunView(APIView):
    """
    Run CustomListQueryTask demo.

    This endpoint reads fixed demo input from:

        DEMO_INPUT_DATA_HOME/custom_list_query/custom_list_query_demo_input.json

    It does not accept user-provided RNA input.
    """

    parser_classes = [JSONParser]

    def post(self, request):
        try:
            try:
                demo_input = load_custom_list_query_demo_input()
            except FileNotFoundError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )
            except (
                CustomListQueryDemoPathError,
                CustomListQueryDemoConfigError,
            ) as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            task_name = str(
                demo_input.get("task_name", "")
            ).strip()

            map_info = str(
                demo_input.get("map_info", "")
            ).strip()

            rnas = demo_input.get("rnas")

            if not task_name:
                return Response(
                    {
                        "success": False,
                        "msg": "Demo input is missing field: task_name.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not map_info:
                return Response(
                    {
                        "success": False,
                        "msg": "Demo input is missing field: map_info.",
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
                prepare_custom_list_query_workspace(task)
            except (
                OSError,
                CustomListQueryTaskInputError,
                CustomListQueryPathError,
            ) as e:
                task.status = CustomListQueryTask.Status.Failed
                task.finish_time = timezone.now()
                task.save(update_fields=["status", "finish_time"])

                return Response(
                    {
                        "success": False,
                        "msg": f"Failed to prepare demo task workspace: {str(e)}",
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
                    "msg": "Demo task submitted successfully.",
                    "data": {
                        "uuid": str(task.uuid),
                        "task_name": task.task_name,
                        "user": task.user,
                        "status": task.get_status_display(),
                        "create_time": timezone.localtime(
                            task.create_time
                        ).strftime("%Y-%m-%d %H:%M:%S"),
                        "map_info": task.map_info,
                        "is_demo": True,
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


class PairedCohortDemoRunView(APIView):
    """
    Run PairedCohortTask demo.

    This endpoint reads fixed demo input from:
        DEMO_INPUT_DATA_HOME/paired_cohort/paired_cohort_demo_input.json

    It creates a new PairedCohortTask, copies demo CSV files into
    task workspace/input, validates input files, and submits the task to Slurm.

    Optional request body:
        {
            "user": "demo_user"
        }
    """

    parser_classes = [JSONParser]

    def post(self, request):
        try:
            try:
                demo_input = load_paired_cohort_demo_input()

            except FileNotFoundError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            except (
                PairedCohortDemoPathError,
                PairedCohortDemoManifestError,
                PairedCohortDemoConfigError,
            ) as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            task_name = str(
                demo_input.get("task_name", "")
            ).strip()

            map_info = str(
                demo_input.get("map_info", "")
            ).strip()

            deg_method = str(
                demo_input.get("deg_method", "")
            ).strip()

            if not task_name:
                return Response(
                    {
                        "success": False,
                        "msg": "Demo input is missing field: task_name.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                validate_task_name_for_filename(task_name)
            except PairedCohortTaskPathError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not map_info:
                return Response(
                    {
                        "success": False,
                        "msg": "Demo input is missing field: map_info.",
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

            if deg_method not in ["limma", "deseq2"]:
                return Response(
                    {
                        "success": False,
                        "msg": (
                            "Invalid demo deg_method. "
                            "Allowed values: limma, deseq2."
                        ),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                cutoff_fields = get_paired_cohort_demo_cutoff_fields(
                    demo_input
                )
            except PairedCohortDemoConfigError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            task = PairedCohortTask.objects.create(
                user=str(
                    request.data.get("user", "")
                ).strip(),
                task_name=task_name,
                status=PairedCohortTask.Status.Pending,
                map_info=map_info,
                deg_method=deg_method,
                **cutoff_fields,
            )

            try:
                prepare_paired_cohort_workspace(task)

                saved_files = copy_paired_cohort_demo_input_files_to_task(
                    task=task,
                    config=demo_input,
                )

                validate_paired_cohort_file_contents(task)

                task.mrna_file = saved_files["mrna_file"]
                task.mirna_file = saved_files["mirna_file"]
                task.lncrna_file = saved_files["lncrna_file"]
                task.circrna_file = saved_files["circrna_file"]
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
                FileNotFoundError,
                PairedCohortTaskPathError,
                PairedCohortDemoPathError,
                PairedCohortDemoConfigError,
            ) as e:
                task.status = PairedCohortTask.Status.Failed
                task.finish_time = timezone.now()
                task.save(update_fields=["status", "finish_time"])

                return Response(
                    {
                        "success": False,
                        "msg": (
                            "Failed to prepare demo task workspace: "
                            f"{str(e)}"
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
                    "msg": "Demo task submitted successfully.",
                    "data": {
                        "uuid": str(task.uuid),
                        "task_type": "PairedCohortTask",
                        "task_name": task.task_name,
                        "user": task.user,
                        "status": task.get_status_display(),
                        "create_time": timezone.localtime(
                            task.create_time
                        ).strftime("%Y-%m-%d %H:%M:%S"),
                        "map_info": task.map_info,
                        "deg_method": task.deg_method,
                        "is_demo": True,
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


class PairedCohortDemoInfoView(APIView):
    """
    Return paired cohort demo manifest and file availability.
    """

    def get(self, request):
        try:
            manifest = load_paired_cohort_demo_manifest()
            demo_dir = get_paired_cohort_demo_dir()

            csv_files = manifest.get("csv_files", {})
            parquet_files = manifest.get("parquet_files", {})

            file_status = {
                "csv_files": {},
                "parquet_files": {},
            }

            for key, filename in csv_files.items():
                file_path = demo_dir / filename
                file_status["csv_files"][key] = {
                    "filename": filename,
                    "exists": file_path.exists() and file_path.is_file(),
                }

            for key, filename in parquet_files.items():
                file_path = demo_dir / filename
                file_status["parquet_files"][key] = {
                    "filename": filename,
                    "exists": file_path.exists() and file_path.is_file(),
                }

            return Response(
                {
                    "workflow_type": "paired_cohort",
                    "task_type": manifest.get("task_type"),
                    "description": manifest.get("description", ""),
                    "valid_rna_types": PAIRED_COHORT_DEMO_VALID_RNA_TYPES,
                    "csv_files": csv_files,
                    "parquet_files": parquet_files,
                    "file_status": file_status,
                },
                status=status.HTTP_200_OK,
            )

        except FileNotFoundError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_404_NOT_FOUND,
            )
        except (
            PairedCohortDemoPathError,
            PairedCohortDemoManifestError,
        ) as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )
        except Exception as e:
            print(traceback.format_exc())
            return Response(
                {"detail": f"Failed to read paired cohort demo info: {e}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )


class PairedCohortDemoSampleMetaView(APIView):
    """
    Return paired cohort demo sample metadata.
    """

    def get(self, request):
        try:
            file_path = get_paired_cohort_demo_meta_file_path()
            validate_demo_file_exists(
                file_path=file_path,
                file_label="Paired cohort demo sample metadata file",
            )

        except FileNotFoundError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_404_NOT_FOUND,
            )
        except (
            PairedCohortDemoPathError,
            PairedCohortDemoManifestError,
        ) as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )

        try:
            df = pd.read_csv(file_path, low_memory=False, encoding="utf-8-sig")
        except Exception as e:
            return Response(
                {"detail": f"Failed to read demo sample metadata file: {e}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        df = df.fillna("")

        return Response(
            {
                "workflow_type": "paired_cohort",
                "file_format": "csv",
                "count": len(df),
                "columns": df.columns.tolist(),
                "results": df.to_dict("records"),
            },
            status=status.HTTP_200_OK,
        )


class PairedCohortDemoExpressionGeneListView(APIView):
    """
    Return gene list from demo expression parquet.
    """

    def get(self, request):
        rna_type = str(
            request.query_params.get("rna_type", "")
        ).strip()

        if not rna_type:
            return Response(
                {"detail": "Missing query parameter: rna_type."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        if rna_type not in PAIRED_COHORT_DEMO_VALID_RNA_TYPES:
            return Response(
                {
                    "detail": (
                        "Invalid rna_type. Allowed values are: "
                        f"{', '.join(PAIRED_COHORT_DEMO_VALID_RNA_TYPES)}."
                    )
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        try:
            file_path = get_paired_cohort_demo_expression_file_path(
                rna_type=rna_type,
            )

            validate_demo_file_exists(
                file_path=file_path,
                file_label=f"Paired cohort demo {rna_type} expression file",
            )

            columns = read_parquet_columns(file_path)

        except FileNotFoundError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_404_NOT_FOUND,
            )
        except (
            PairedCohortDemoError,
            PairedCohortDemoPathError,
            PairedCohortDemoManifestError,
        ) as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )
        except Exception as e:
            print(traceback.format_exc())
            return Response(
                {"detail": f"Failed to read demo expression schema: {e}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        if not columns:
            return Response(
                {"detail": "Demo expression file has no columns."},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        sample_column = columns[0]
        genes = columns[1:]

        return Response(
            {
                "workflow_type": "paired_cohort",
                "rna_type": rna_type,
                "file_format": DEFAULT_DEMO_EXPRESSION_FILE_FORMAT,
                "sample_column": sample_column,
                "count": len(genes),
                "genes": genes,
            },
            status=status.HTTP_200_OK,
        )


class PairedCohortDemoExpressionDataView(APIView):
    """
    Return selected gene expression data from demo expression parquet.

    Body:
        {
            "rna_type": "mRNA",
            "genes": ["GENE1", "GENE2"]
        }
    """

    def post(self, request):
        rna_type = str(
            request.data.get("rna_type", "")
        ).strip()

        genes = request.data.get("genes", [])

        if not rna_type:
            return Response(
                {"detail": "Missing field: rna_type."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        if rna_type not in PAIRED_COHORT_DEMO_VALID_RNA_TYPES:
            return Response(
                {
                    "detail": (
                        "Invalid rna_type. Allowed values are: "
                        f"{', '.join(PAIRED_COHORT_DEMO_VALID_RNA_TYPES)}."
                    )
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        if not isinstance(genes, list):
            return Response(
                {"detail": "Field 'genes' must be a list."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        genes = [
            str(gene).strip()
            for gene in genes
            if str(gene).strip()
        ]

        if not genes:
            return Response(
                {"detail": "At least one gene must be selected."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        if len(genes) > PAIRED_COHORT_DEMO_MAX_SELECTED_GENES:
            return Response(
                {
                    "detail": (
                        "At most "
                        f"{PAIRED_COHORT_DEMO_MAX_SELECTED_GENES} "
                        "genes can be selected."
                    )
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        try:
            file_path = get_paired_cohort_demo_expression_file_path(
                rna_type=rna_type,
            )

            validate_demo_file_exists(
                file_path=file_path,
                file_label=f"Paired cohort demo {rna_type} expression file",
            )

            all_columns = read_parquet_columns(file_path)

        except FileNotFoundError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_404_NOT_FOUND,
            )
        except (
            PairedCohortDemoError,
            PairedCohortDemoPathError,
            PairedCohortDemoManifestError,
        ) as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )
        except Exception as e:
            print(traceback.format_exc())
            return Response(
                {"detail": f"Failed to read demo expression schema: {e}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        if not all_columns:
            return Response(
                {"detail": "Demo expression file has no columns."},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        sample_column = all_columns[0]
        available_genes = set(all_columns[1:])

        missing_genes = [
            gene for gene in genes
            if gene not in available_genes
        ]

        if missing_genes:
            return Response(
                {
                    "detail": "Some genes are not found in the demo expression file.",
                    "missing_genes": missing_genes,
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        usecols = [sample_column] + genes

        try:
            df = pd.read_parquet(
                file_path,
                columns=usecols,
                engine="pyarrow",
            )
        except Exception as e:
            return Response(
                {"detail": f"Failed to read demo expression data: {e}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        df = df.fillna("")

        return Response(
            {
                "workflow_type": "paired_cohort",
                "rna_type": rna_type,
                "file_format": DEFAULT_DEMO_EXPRESSION_FILE_FORMAT,
                "count": len(df),
                "columns": df.columns.tolist(),
                "results": df.to_dict("records"),
            },
            status=status.HTTP_200_OK,
        )


class PairedCohortDemoDataDownloadView(APIView):
    """
    Download paired cohort demo input data zip.

    The zip file must already exist under paired cohort demo input directory:

        paired_cohort_demo_data.zip

    GET /api/analysis/paired_cohort_demo/download_data/
    """

    def get(self, request):
        try:
            archive_path = validate_paired_cohort_demo_data_archive()

            return FileResponse(
                open(archive_path, "rb"),
                as_attachment=True,
                filename=PAIRED_COHORT_DEMO_DATA_ARCHIVE_NAME,
                content_type="application/zip",
            )

        except FileNotFoundError as e:
            return Response(
                {
                    "success": False,
                    "detail": str(e),
                },
                status=status.HTTP_404_NOT_FOUND,
            )

        except PairedCohortDemoPathError as e:
            return Response(
                {
                    "success": False,
                    "detail": str(e),
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {
                    "success": False,
                    "detail": (
                        "Failed to download paired cohort demo data archive: "
                        f"{str(e)}"
                    ),
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )
