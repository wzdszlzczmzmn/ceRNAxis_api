import traceback

from django.conf import settings

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.utils.workflow_detail_utils.workflow_log2fc_background_utils import (
    PAIRED_COHORT_VALID_BACKGROUND_TYPES,
    HYBRID_REFERENCE_VALID_BACKGROUND_TYPES,
    WORKFLOW_LOG2FC_X_COL,
    WORKFLOW_LOG2FC_Y_COL,
    WorkflowLog2FCBackgroundInputError,
    read_log2fc_background_file_by_path,
    get_workflow_available_background_types,
    build_log2fc_correlation_response_data_from_dataframe,
)
from database.utils.dataset_annotation_utils.path_utils import (
    DatasetAnnotationInputError,
    DatasetAnnotationPathError,
    get_dataset_query_name,
    resolve_dataset_annotation_dir,
    resolve_tcga_annotation_dir_name,
    resolve_timedb_annotation_dir_name,
    get_dataset_annotation_log2fc_background_file_path,
)


class BaseDatasetAnnotationLog2FCCorrelationView(APIView):
    """
    Return dataset-level annotation log2FC background correlation plot data.

    Query params:
        dataset or dataset_name
        type: optional interaction type

    Input filename:
        {annotation_file_prefix}_ceRNA_background.csv

    Plot mapping:
        x = ceRNA_log2FC
        y = miRNA_log2FC
    """

    source = None
    network_source_task_type = None

    annotation_root_setting_name = None
    annotation_dir_name_resolver = None

    valid_background_types = []
    x_col = WORKFLOW_LOG2FC_X_COL
    y_col = WORKFLOW_LOG2FC_Y_COL

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

    def resolve_annotation_context(self, dataset_name: str) -> dict:
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

        background_file = get_dataset_annotation_log2fc_background_file_path(
            annotation_dir=annotation_dir,
            file_prefix=file_prefix,
        )

        return {
            "dataset_name": dataset_name,
            "annotation_dir_name": annotation_dir_name,
            "annotation_file_prefix": file_prefix,
            "annotation_dir": annotation_dir,
            "background_file": background_file,
        }

    def get(self, request):
        try:
            if not self.source:
                raise RuntimeError("Missing source.")

            if not self.network_source_task_type:
                raise RuntimeError("Missing network_source_task_type.")

            if not self.valid_background_types:
                raise RuntimeError("Missing valid_background_types.")

            dataset_name = get_dataset_query_name(request)

            requested_type = request.query_params.get("type", None)

            type_value = (
                str(requested_type).strip()
                if requested_type is not None
                else ""
            )

            try:
                context = self.resolve_annotation_context(dataset_name)

                background_file, df = read_log2fc_background_file_by_path(
                    context["background_file"]
                )

            except FileNotFoundError as e:
                return Response(
                    {
                        "detail": str(e),
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            except (
                DatasetAnnotationInputError,
                DatasetAnnotationPathError,
                WorkflowLog2FCBackgroundInputError,
            ) as e:
                return Response(
                    {
                        "detail": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            available_types = get_workflow_available_background_types(
                df=df,
                valid_types=self.valid_background_types,
            )

            if not available_types:
                return Response(
                    {
                        "detail": (
                            "No supported background interaction type found "
                            "in ceRNA background file."
                        )
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not type_value:
                type_value = available_types[0]

            if type_value not in available_types:
                return Response(
                    {
                        "detail": (
                            "Invalid type for this dataset annotation. "
                            "Allowed values are: "
                            f"{', '.join(available_types)}."
                        ),
                        "available_types": available_types,
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            base_response = {
                "success": True,
                "source": self.source,
                "dataset_name": context["dataset_name"],
                "annotation_dir_name": context["annotation_dir_name"],
                "annotation_file_prefix": context[
                    "annotation_file_prefix"
                ],
                "network_source_task_type": self.network_source_task_type,
            }

            result = build_log2fc_correlation_response_data_from_dataframe(
                df=df,
                type_value=type_value,
                background_file_name=background_file.name,
                available_types=available_types,
                base_response=base_response,
                x_col=self.x_col,
                y_col=self.y_col,
            )

            return Response(
                result,
                status=status.HTTP_200_OK,
            )

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {
                    "detail": f"Server error: {str(e)}",
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )


class TCGADatasetAnnotationLog2FCCorrelationView(
    BaseDatasetAnnotationLog2FCCorrelationView
):
    """
    TCGA dataset annotation log2FC background correlation.

    Source semantics:
        Paired Cohort annotation output.

    Valid interaction types:
        miRNA-mRNA, miRNA-lncRNA, miRNA-circRNA
    """

    source = "TCGA"
    network_source_task_type = "PairedCohortTask"

    annotation_root_setting_name = "TCGA_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_tcga_annotation_dir_name
    )

    valid_background_types = PAIRED_COHORT_VALID_BACKGROUND_TYPES


class TIMEDBDatasetAnnotationLog2FCCorrelationView(
    BaseDatasetAnnotationLog2FCCorrelationView
):
    """
    TIMEDB dataset annotation log2FC background correlation.

    Source semantics:
        Hybrid Reference annotation output.

    Valid interaction types:
        miRNA-mRNA, miRNA-lncRNA
    """

    source = "TIMEDB"
    network_source_task_type = "HybridReferenceTask"

    annotation_root_setting_name = "TIMEDB_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_timedb_annotation_dir_name
    )

    valid_background_types = HYBRID_REFERENCE_VALID_BACKGROUND_TYPES
