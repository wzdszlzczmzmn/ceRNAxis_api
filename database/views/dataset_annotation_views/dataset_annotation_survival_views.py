import traceback

import pandas as pd

from django.conf import settings
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.utils.paired_cohort_task_utils import (
    PAIRED_COHORT_SURVIVAL_GROUPS,
)
from analysis.utils.hybrid_reference_task_utils import (
    HYBRID_REFERENCE_SURVIVAL_GROUPS,
)
from analysis.utils.workflow_detail_utils.survival_km_utils import (
    SurvivalKMInputError,
    validate_survival_dataframe_columns,
    build_survival_km_data_from_dataframe_common,
)

from database.utils.dataset_annotation_utils.path_utils import (
    DatasetAnnotationInputError,
    DatasetAnnotationPathError,
    get_dataset_query_name,
    resolve_dataset_annotation_dir,
    resolve_tcga_annotation_dir_name,
    resolve_timedb_annotation_dir_name,
    get_dataset_annotation_survival_file_path,
)


class BaseDatasetAnnotationSurvivalKMDataView(APIView):
    """
    Return Kaplan-Meier survival curve data for dataset annotation.

    Query params:
        dataset or dataset_name

    Input filename:
        {annotation_file_prefix}_survival_analysis.csv
    """

    source = None
    network_source_task_type = None

    annotation_root_setting_name = None
    annotation_dir_name_resolver = None

    title = "ceRNA axis-based survival analysis"
    valid_groups = []

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

    def resolve_annotation_context(self, request) -> dict:
        dataset_name = get_dataset_query_name(request)

        annotation_dir_name = self.annotation_dir_name_resolver(
            dataset_name
        )

        annotation_dir = resolve_dataset_annotation_dir(
            annotation_root_dir=self.get_annotation_root_dir(),
            annotation_dir_name=annotation_dir_name,
        )

        file_prefix = self.get_annotation_file_prefix(
            dataset_name=dataset_name,
            annotation_dir_name=annotation_dir_name,
        )

        survival_file = get_dataset_annotation_survival_file_path(
            annotation_dir=annotation_dir,
            file_prefix=file_prefix,
        )

        return {
            "dataset_name": dataset_name,
            "annotation_dir_name": annotation_dir_name,
            "annotation_file_prefix": file_prefix,
            "annotation_dir": annotation_dir,
            "survival_file": survival_file,
        }

    def read_survival_file(self, survival_file):
        if not survival_file.exists() or not survival_file.is_file():
            raise FileNotFoundError(
                f"Dataset annotation survival analysis file not found: "
                f"{survival_file.name}"
            )

        try:
            df = pd.read_csv(survival_file)
        except Exception as e:
            raise SurvivalKMInputError(
                f"Failed to read dataset annotation survival analysis file: "
                f"{str(e)}"
            )

        validate_survival_dataframe_columns(df=df)

        return df

    def get_base_response(self, context: dict) -> dict:
        return {
            "success": True,
            "source": self.source,
            "dataset_name": context["dataset_name"],
            "annotation_dir_name": context["annotation_dir_name"],
            "annotation_file_prefix": context["annotation_file_prefix"],
            "network_source_task_type": self.network_source_task_type,
        }

    def get(self, request):
        try:
            if not self.source:
                raise RuntimeError("Missing source.")

            if not self.network_source_task_type:
                raise RuntimeError("Missing network_source_task_type.")

            if not self.annotation_dir_name_resolver:
                raise RuntimeError("Missing annotation_dir_name_resolver.")

            if not self.valid_groups:
                raise RuntimeError("Missing valid_groups.")

            context = self.resolve_annotation_context(request)

            try:
                df = self.read_survival_file(
                    context["survival_file"]
                )

                result = build_survival_km_data_from_dataframe_common(
                    survival_file_name=context["survival_file"].name,
                    df=df,
                    title=self.title,
                    base_response=self.get_base_response(context),
                    valid_groups=self.valid_groups,
                )

            except FileNotFoundError as e:
                return Response(
                    {"success": False, "detail": str(e)},
                    status=status.HTTP_404_NOT_FOUND,
                )

            except SurvivalKMInputError as e:
                return Response(
                    {"success": False, "detail": str(e)},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            return Response(
                result,
                status=status.HTTP_200_OK,
            )

        except DatasetAnnotationInputError as e:
            return Response(
                {"success": False, "detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )

        except DatasetAnnotationPathError as e:
            return Response(
                {"success": False, "detail": str(e)},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
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


class TCGADatasetAnnotationSurvivalKMDataView(
    BaseDatasetAnnotationSurvivalKMDataView
):
    """
    TCGA dataset annotation survival KM data.

    Source semantics:
        Paired Cohort annotation output.

    Input filename:
        {annotation_file_prefix}_survival_analysis.csv
    """

    source = "TCGA"
    network_source_task_type = "PairedCohortTask"

    annotation_root_setting_name = "TCGA_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_tcga_annotation_dir_name
    )

    title = "ceRNA axis-based survival analysis"
    valid_groups = PAIRED_COHORT_SURVIVAL_GROUPS


class TIMEDBDatasetAnnotationSurvivalKMDataView(
    BaseDatasetAnnotationSurvivalKMDataView
):
    """
    TIMEDB dataset annotation survival KM data.

    Source semantics:
        Hybrid Reference annotation output.

    Input filename:
        {annotation_file_prefix}_survival_analysis.csv
    """

    source = "TIMEDB"
    network_source_task_type = "HybridReferenceTask"

    annotation_root_setting_name = "TIMEDB_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_timedb_annotation_dir_name
    )

    title = "TCGA-based ceRNA axis survival analysis"
    valid_groups = HYBRID_REFERENCE_SURVIVAL_GROUPS
