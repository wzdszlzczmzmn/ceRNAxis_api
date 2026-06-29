import traceback

from django.conf import settings

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.utils.workflow_detail_utils.workflow_axis_final_utils import (
    WorkflowAxisFinalInputError,
    build_axis_final_response_from_dataframe,
    read_axis_final_file_by_path,
    PAIRED_COHORT_AXIS_FINAL_COLUMNS,
    HYBRID_REFERENCE_AXIS_FINAL_COLUMNS,
    WORKFLOW_AXIS_FINAL_NUMERIC_COLUMNS, PAIRED_COHORT_AXIS_FINAL_REQUIRED_COLUMNS,
    HYBRID_REFERENCE_AXIS_FINAL_REQUIRED_COLUMNS,
)
from database.utils.dataset_annotation_utils.path_utils import (
    DatasetAnnotationInputError,
    DatasetAnnotationPathError,
    get_dataset_query_name,
    resolve_dataset_annotation_dir,
    get_dataset_annotation_axis_final_file_path,
    resolve_tcga_annotation_dir_name,
    resolve_timedb_annotation_dir_name,
)


class BaseDatasetAnnotationAxisFinalDataView(APIView):
    source = None
    network_source_task_type = None
    annotation_root_setting_name = None
    annotation_dir_name_resolver = None

    axis_final_columns = []
    axis_final_required_columns = []
    axis_final_numeric_columns = WORKFLOW_AXIS_FINAL_NUMERIC_COLUMNS

    def get_annotation_root_dir(self):
        if not self.annotation_root_setting_name:
            raise DatasetAnnotationPathError(
                "Annotation root setting name is not configured."
            )

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

    def get_annotation_file_prefix(
        self,
        *,
        dataset_name: str,
        annotation_dir_name: str,
    ) -> str:
        return annotation_dir_name

    def get_axis_final_file_info(self, dataset_name: str) -> dict:
        if self.annotation_dir_name_resolver is None:
            raise DatasetAnnotationPathError(
                "Annotation directory resolver is not configured."
            )

        annotation_dir_name = self.annotation_dir_name_resolver(dataset_name)

        annotation_dir = resolve_dataset_annotation_dir(
            annotation_root_dir=self.get_annotation_root_dir(),
            annotation_dir_name=annotation_dir_name,
        )

        file_prefix = self.get_annotation_file_prefix(
            dataset_name=dataset_name,
            annotation_dir_name=annotation_dir_name,
        )

        axis_final_file = get_dataset_annotation_axis_final_file_path(
            annotation_dir=annotation_dir,
            file_prefix=file_prefix,
        )

        return {
            "dataset_name": dataset_name,
            "annotation_dir_name": annotation_dir_name,
            "annotation_file_prefix": file_prefix,
            "axis_final_file": axis_final_file,
        }

    def get(self, request):
        try:
            if not self.source:
                raise RuntimeError("Missing source.")

            if not self.network_source_task_type:
                raise RuntimeError("Missing network_source_task_type.")

            if not self.axis_final_columns:
                raise RuntimeError("Missing axis_final_columns.")

            dataset_name = get_dataset_query_name(request)

            required_columns = (
                self.axis_final_required_columns
                or self.axis_final_columns
            )

            try:
                file_info = self.get_axis_final_file_info(dataset_name)

                axis_file, df = read_axis_final_file_by_path(
                    file_path=file_info["axis_final_file"],
                    required_columns=required_columns,
                )

            except FileNotFoundError as e:
                return Response(
                    {
                        "success": False,
                        "detail": str(e),
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            except (
                DatasetAnnotationInputError,
                DatasetAnnotationPathError,
                WorkflowAxisFinalInputError,
            ) as e:
                return Response(
                    {
                        "success": False,
                        "detail": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            base_response = {
                "success": True,
                "source": self.source,
                "dataset_name": file_info["dataset_name"],
                "annotation_dir_name": file_info["annotation_dir_name"],
                "annotation_file_prefix": file_info[
                    "annotation_file_prefix"
                ],
                "network_source_task_type": self.network_source_task_type,
            }

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

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {
                    "success": False,
                    "detail": f"Server error: {str(e)}",
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )


class TCGADatasetAnnotationAxisFinalDataView(
    BaseDatasetAnnotationAxisFinalDataView
):
    """
    TCGA dataset annotation axis final data.

    Input:
        ?dataset=TCGA_ACC_mRNA

    Resolution:
        TCGA_ACC_mRNA -> TCGA_ACC

    Source semantics:
        Paired Cohort annotation output.
    """

    source = "TCGA"
    network_source_task_type = "PairedCohortTask"

    annotation_root_setting_name = "TCGA_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_tcga_annotation_dir_name
    )

    axis_final_columns = PAIRED_COHORT_AXIS_FINAL_COLUMNS
    axis_final_required_columns = PAIRED_COHORT_AXIS_FINAL_REQUIRED_COLUMNS


class TIMEDBDatasetAnnotationAxisFinalDataView(
    BaseDatasetAnnotationAxisFinalDataView
):
    """
    TIMEDB dataset annotation axis final data.

    Source semantics:
        Hybrid Reference annotation output.
    """

    source = "TIMEDB"
    network_source_task_type = "HybridReferenceTask"

    annotation_root_setting_name = "TIMEDB_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_timedb_annotation_dir_name
    )

    axis_final_columns = HYBRID_REFERENCE_AXIS_FINAL_COLUMNS
    axis_final_required_columns = HYBRID_REFERENCE_AXIS_FINAL_REQUIRED_COLUMNS
