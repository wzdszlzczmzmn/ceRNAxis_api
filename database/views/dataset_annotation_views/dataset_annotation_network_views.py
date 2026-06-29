import traceback

from django.conf import settings

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.utils.workflow_detail_utils.workflow_output_network_builder import (
    WorkflowOutputNetworkBuilderMixin,
)
from database.utils.dataset_annotation_utils.path_utils import (
    DatasetAnnotationInputError,
    DatasetAnnotationPathError,
    get_dataset_query_name,
    resolve_dataset_annotation_dir,
    get_dataset_annotation_file_path,
    resolve_tcga_annotation_dir_name,
    resolve_timedb_annotation_dir_name,
)


DATASET_ANNOTATION_CERNA_AXIS_SUFFIX = "_ceRNA_axis.csv"
DATASET_ANNOTATION_IMMUNE_AXIS_SUFFIX = "_map_immune_axis.csv"


class BaseDatasetAnnotationNetworkView(
    WorkflowOutputNetworkBuilderMixin,
    APIView,
):
    source = None
    network_source_task_type = None
    annotation_root_setting_name = None
    annotation_dir_name_resolver = None

    cerna_node_source = "dataset_annotation_cerna_axis"
    cerna_edge_source = "dataset_annotation_cerna_axis"
    immune_node_source = "dataset_annotation_immune_axis"
    immune_edge_source = "dataset_annotation_immune_axis"

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
        默认认为 annotation 文件名前缀与 annotation 目录名一致。

        Example:
            TCGA_ACC/
                TCGA_ACC_ceRNA_axis.csv
                TCGA_ACC_map_immune_axis.csv
        """
        return annotation_dir_name

    def get_network_files(self, dataset_name: str) -> dict:
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

        ceRNA_axis_file = get_dataset_annotation_file_path(
            annotation_dir=annotation_dir,
            file_prefix=file_prefix,
            filename_suffix=DATASET_ANNOTATION_CERNA_AXIS_SUFFIX,
        )

        immune_axis_file = get_dataset_annotation_file_path(
            annotation_dir=annotation_dir,
            file_prefix=file_prefix,
            filename_suffix=DATASET_ANNOTATION_IMMUNE_AXIS_SUFFIX,
        )

        return {
            "dataset_name": dataset_name,
            "annotation_dir_name": annotation_dir_name,
            "annotation_file_prefix": file_prefix,
            "annotation_dir": annotation_dir,
            "ceRNA_axis_file": ceRNA_axis_file,
            "immune_axis_file": immune_axis_file,
        }

    def get(self, request, *args, **kwargs):
        try:
            dataset_name = get_dataset_query_name(request)

            try:
                file_info = self.get_network_files(dataset_name)

                result = self.build_network_from_files(
                    ceRNA_axis_file=file_info["ceRNA_axis_file"],
                    immune_axis_file=file_info["immune_axis_file"],
                )

            except (
                DatasetAnnotationInputError,
                DatasetAnnotationPathError,
                ValueError,
            ) as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            result.update(
                {
                    "success": True,
                    "source": self.source,
                    "dataset_name": file_info["dataset_name"],
                    "annotation_dir_name": file_info["annotation_dir_name"],
                    "annotation_file_prefix": file_info[
                        "annotation_file_prefix"
                    ],
                    "network_source_task_type": (
                        self.network_source_task_type
                    ),
                }
            )

            return Response(result, status=status.HTTP_200_OK)

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {
                    "success": False,
                    "msg": f"Server error: {str(e)}",
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )


class TCGADatasetAnnotationNetworkView(
    BaseDatasetAnnotationNetworkView
):
    """
    TCGA dataset annotation network.

    Input:
        ?dataset=TCGA_ACC_mRNA

    Resolution:
        TCGA_ACC_mRNA -> TCGA_ACC

    File expectation:
        TCGA_DATASET_ANNOTATIONS_DIR/
            TCGA_ACC/
                TCGA_ACC_ceRNA_axis.csv
                TCGA_ACC_map_immune_axis.csv

    Annotation source:
        Paired Cohort workflow output.
    """

    source = "TCGA"
    network_source_task_type = "PairedCohortTask"
    annotation_root_setting_name = "TCGA_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_tcga_annotation_dir_name
    )

    cerna_node_source = "paired_cohort_cerna_axis"
    cerna_edge_source = "paired_cohort_cerna_axis"
    immune_node_source = "paired_cohort_immune_axis"
    immune_edge_source = "paired_cohort_immune_axis"


class TIMEDBDatasetAnnotationNetworkView(
    BaseDatasetAnnotationNetworkView
):
    """
    TIMEDB dataset annotation network.

    Input:
        ?dataset=GSE19750_mRNA or ?dataset=GSE19750

    Annotation source:
        Hybrid Reference workflow output.
    """

    source = "TIMEDB"
    network_source_task_type = "HybridReferenceTask"
    annotation_root_setting_name = "TIMEDB_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_timedb_annotation_dir_name
    )

    cerna_node_source = "hybrid_reference_cerna_axis"
    cerna_edge_source = "hybrid_reference_cerna_axis"
    immune_node_source = "hybrid_reference_immune_axis"
    immune_edge_source = "hybrid_reference_immune_axis"
