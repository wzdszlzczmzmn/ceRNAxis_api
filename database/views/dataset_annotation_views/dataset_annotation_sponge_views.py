import traceback

from django.conf import settings

from rest_framework import status
from rest_framework.response import Response
from rest_framework.views import APIView

from analysis.utils.workflow_detail_utils.workflow_sponge_utils import (
    WORKFLOW_SPONGE_COLUMNS,
    WORKFLOW_SPONGE_NUMERIC_COLUMNS,
    WORKFLOW_SPONGE_REQUIRED_COLUMNS,
    WorkflowSpongeInputError,
    build_sponge_response_from_dataframe,
    read_sponge_file_by_path,
)

from database.utils.dataset_annotation_utils.path_utils import (
    DatasetAnnotationInputError,
    DatasetAnnotationPathError,
    get_dataset_annotation_sponge_result_file_path,
    get_dataset_query_name,
    get_timedb_group_type_query,
    resolve_dataset_annotation_dir,
    resolve_tcga_annotation_dir_name,
    resolve_timedb_group_annotation_dir_name,
    resolve_timedb_group_annotation_file_prefix,
)


class TCGADatasetAnnotationSpongeDataView(APIView):
    """
    Return Sponge results for a TCGA dataset annotation.

    Query parameters:
        dataset:
            Dataset name, for example TCGA_ACC_mRNA.

    Resolution:
        TCGA_ACC_mRNA -> TCGA_ACC

    Expected file:
        {TCGA_DATASET_ANNOTATIONS_DIR}/
            TCGA_ACC/
                TCGA_ACC_sponge_result.csv

    Source semantics:
        Paired Cohort annotation output.
    """

    source = "TCGA"
    network_source_task_type = "PairedCohortTask"
    annotation_root_setting_name = "TCGA_DATASET_ANNOTATIONS_DIR"

    sponge_columns = WORKFLOW_SPONGE_COLUMNS
    sponge_required_columns = WORKFLOW_SPONGE_REQUIRED_COLUMNS
    sponge_numeric_columns = WORKFLOW_SPONGE_NUMERIC_COLUMNS

    def get_annotation_root_dir(self):
        annotation_root_dir = getattr(
            settings,
            self.annotation_root_setting_name,
            None,
        )

        if not annotation_root_dir:
            raise DatasetAnnotationPathError(
                f"{self.annotation_root_setting_name} is not configured."
            )

        return annotation_root_dir

    def get(self, request):
        try:
            dataset_name = get_dataset_query_name(request)

            annotation_dir_name = resolve_tcga_annotation_dir_name(
                dataset_name
            )

            annotation_dir = resolve_dataset_annotation_dir(
                annotation_root_dir=self.get_annotation_root_dir(),
                annotation_dir_name=annotation_dir_name,
            )

            file_prefix = annotation_dir_name

            sponge_file_path = (
                get_dataset_annotation_sponge_result_file_path(
                    annotation_dir=annotation_dir,
                    file_prefix=file_prefix,
                )
            )

            sponge_file, df = read_sponge_file_by_path(
                file_path=sponge_file_path,
                required_columns=self.sponge_required_columns,
            )

            response_data = build_sponge_response_from_dataframe(
                df=df,
                sponge_file_name=sponge_file.name,
                columns=self.sponge_columns,
                required_columns=self.sponge_required_columns,
                numeric_columns=self.sponge_numeric_columns,
                base_response={
                    "success": True,
                    "source": self.source,
                    "dataset_name": dataset_name,
                    "network_source_task_type": (
                        self.network_source_task_type
                    ),
                },
            )

            return Response(
                response_data,
                status=status.HTTP_200_OK,
            )

        except DatasetAnnotationInputError as exc:
            return Response(
                {
                    "success": False,
                    "detail": str(exc),
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        except FileNotFoundError as exc:
            return Response(
                {
                    "success": False,
                    "detail": str(exc),
                },
                status=status.HTTP_404_NOT_FOUND,
            )

        except WorkflowSpongeInputError as exc:
            return Response(
                {
                    "success": False,
                    "detail": str(exc),
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        except DatasetAnnotationPathError as exc:
            return Response(
                {
                    "success": False,
                    "detail": str(exc),
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        except Exception as exc:
            print(traceback.format_exc())

            return Response(
                {
                    "success": False,
                    "detail": f"Server error: {str(exc)}",
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )


class TIMEDBDatasetAnnotationSpongeDataView(APIView):
    """
    Return Sponge results for one TIMEDB dataset annotation group.

    Query parameters:
        dataset:
            Dataset name, for example GSE20194.

        group_type:
            other | grade | stage

        group_by:
            Current frontend group-by value. It is returned as response
            context, but file resolution is controlled by group_type.

    Expected files:
        group_type=other:
            {TIMEDB_DATASET_ANNOTATIONS_DIR}/
                GSE20194/
                    GSE20194_sponge_result.csv

        group_type=grade:
            {TIMEDB_DATASET_ANNOTATIONS_DIR}/
                GSE20194_grade/
                    GSE20194_sponge_result.csv

        group_type=stage:
            {TIMEDB_DATASET_ANNOTATIONS_DIR}/
                GSE20194_stage/
                    GSE20194_sponge_result.csv

    Source semantics:
        Hybrid Reference annotation output.
    """

    source = "TIMEDB"
    network_source_task_type = "HybridReferenceTask"
    annotation_root_setting_name = "TIMEDB_DATASET_ANNOTATIONS_DIR"

    sponge_columns = WORKFLOW_SPONGE_COLUMNS
    sponge_required_columns = WORKFLOW_SPONGE_REQUIRED_COLUMNS
    sponge_numeric_columns = WORKFLOW_SPONGE_NUMERIC_COLUMNS

    def get_annotation_root_dir(self):
        annotation_root_dir = getattr(
            settings,
            self.annotation_root_setting_name,
            None,
        )

        if not annotation_root_dir:
            raise DatasetAnnotationPathError(
                f"{self.annotation_root_setting_name} is not configured."
            )

        return annotation_root_dir

    def get_group_by(self, request):
        group_by = request.query_params.get("group_by")

        if group_by is None:
            return None

        group_by = str(group_by).strip()

        return group_by or None

    def get(self, request):
        try:
            dataset_name = get_dataset_query_name(request)
            group_type = get_timedb_group_type_query(request)
            group_by = self.get_group_by(request)

            annotation_dir_name = (
                resolve_timedb_group_annotation_dir_name(
                    dataset_name=dataset_name,
                    group_type=group_type,
                )
            )

            annotation_dir = resolve_dataset_annotation_dir(
                annotation_root_dir=self.get_annotation_root_dir(),
                annotation_dir_name=annotation_dir_name,
            )

            file_prefix = (
                resolve_timedb_group_annotation_file_prefix(
                    dataset_name=dataset_name,
                    group_type=group_type,
                )
            )

            sponge_file_path = (
                get_dataset_annotation_sponge_result_file_path(
                    annotation_dir=annotation_dir,
                    file_prefix=file_prefix,
                )
            )

            sponge_file, df = read_sponge_file_by_path(
                file_path=sponge_file_path,
                required_columns=self.sponge_required_columns,
            )

            response_data = build_sponge_response_from_dataframe(
                df=df,
                sponge_file_name=sponge_file.name,
                columns=self.sponge_columns,
                required_columns=self.sponge_required_columns,
                numeric_columns=self.sponge_numeric_columns,
                base_response={
                    "success": True,
                    "source": self.source,
                    "dataset_name": dataset_name,
                    "group_by": group_by,
                    "group_type": group_type,
                    "network_source_task_type": (
                        self.network_source_task_type
                    ),
                },
            )

            return Response(
                response_data,
                status=status.HTTP_200_OK,
            )

        except DatasetAnnotationInputError as exc:
            return Response(
                {
                    "success": False,
                    "detail": str(exc),
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        except FileNotFoundError as exc:
            return Response(
                {
                    "success": False,
                    "detail": str(exc),
                },
                status=status.HTTP_404_NOT_FOUND,
            )

        except WorkflowSpongeInputError as exc:
            return Response(
                {
                    "success": False,
                    "detail": str(exc),
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        except DatasetAnnotationPathError as exc:
            return Response(
                {
                    "success": False,
                    "detail": str(exc),
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        except Exception as exc:
            print(traceback.format_exc())

            return Response(
                {
                    "success": False,
                    "detail": f"Server error: {str(exc)}",
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )
