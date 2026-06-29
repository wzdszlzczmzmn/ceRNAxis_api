import traceback

from django.conf import settings

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.utils.workflow_detail_utils.workflow_cmap_utils import (
    WORKFLOW_CMAP_REQUIRED_COLUMNS,
    HYBRID_REFERENCE_CMAP_REQUIRED_COLUMNS,
    WorkflowCMapInputError,
    read_cmap_file_by_path,
    build_cmap_response_from_dataframe,
)
from database.utils.dataset_annotation_utils.path_utils import (
    DatasetAnnotationInputError,
    DatasetAnnotationPathError,
    get_dataset_query_name,
    resolve_dataset_annotation_dir,
    get_dataset_annotation_cmap_file_path,
    resolve_tcga_annotation_dir_name,
    resolve_timedb_annotation_dir_name,
)


class BaseDatasetAnnotationCMapResultView(APIView):
    """
    Return dataset-level annotation CMap result table.

    Query params:
        dataset or dataset_name

    Filename:
        {annotation_file_prefix}_CMap.csv
    """

    source = None
    network_source_task_type = None

    annotation_root_setting_name = None
    annotation_dir_name_resolver = None

    required_columns = WORKFLOW_CMAP_REQUIRED_COLUMNS

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
        """
        Default:
            TCGA_ACC/
                TCGA_ACC_CMap.csv
        """
        return annotation_dir_name

    def get_cmap_file_info(self, dataset_name: str) -> dict:
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

        cmap_file = get_dataset_annotation_cmap_file_path(
            annotation_dir=annotation_dir,
            file_prefix=file_prefix,
        )

        return {
            "dataset_name": dataset_name,
            "annotation_dir_name": annotation_dir_name,
            "annotation_file_prefix": file_prefix,
            "cmap_file": cmap_file,
        }

    def get(self, request):
        try:
            if not self.source:
                raise RuntimeError("Missing source.")

            if not self.network_source_task_type:
                raise RuntimeError("Missing network_source_task_type.")

            dataset_name = get_dataset_query_name(request)

            try:
                file_info = self.get_cmap_file_info(dataset_name)

                cmap_file, df = read_cmap_file_by_path(
                    file_path=file_info["cmap_file"],
                    required_columns=self.required_columns,
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
                WorkflowCMapInputError,
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

            response_data = build_cmap_response_from_dataframe(
                df=df,
                cmap_file_name=cmap_file.name,
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


class TCGADatasetAnnotationCMapResultView(
    BaseDatasetAnnotationCMapResultView
):
    """
    TCGA dataset annotation CMap.

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

    # Paired Cohort CMap 当前不强制列结构。
    required_columns = WORKFLOW_CMAP_REQUIRED_COLUMNS


class TIMEDBDatasetAnnotationCMapResultView(
    BaseDatasetAnnotationCMapResultView
):
    """
    TIMEDB dataset annotation CMap.

    Source semantics:
        Hybrid Reference annotation output.
    """

    source = "TIMEDB"
    network_source_task_type = "HybridReferenceTask"

    annotation_root_setting_name = "TIMEDB_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_timedb_annotation_dir_name
    )

    required_columns = HYBRID_REFERENCE_CMAP_REQUIRED_COLUMNS
