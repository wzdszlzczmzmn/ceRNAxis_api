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
    SCST_HYBRID_REFERENCE_SURVIVAL_GROUPS,
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
    get_dataset_annotation_survival_file_path, resolve_timedb_group_annotation_file_prefix,
    resolve_timedb_group_annotation_dir_name, get_timedb_group_type_query,
)

from database.utils.dataset_annotation_utils.scst_path_utils import (
    get_scst_data_type_query,
    get_scst_group_by_query,
    get_scst_group_value_query,
    validate_scst_dataset_group_by,
    resolve_scst_dataset_group_annotation_dir,
    get_scst_dataset_survival_file_path,
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

    # TCGA/TIMEDB use the generic annotation-dir resolver.
    # SC/ST resolves its directory from data_type + group_by.
    requires_annotation_dir_name_resolver = True

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

    def get_group_type(self, request):
        return None

    def get_group_by(self, request):
        group_by = request.query_params.get("group_by")

        if group_by is None:
            return None

        group_by = str(group_by).strip()

        return group_by or None

    def get_annotation_dir_name(
            self,
            *,
            dataset_name: str,
            group_type: str | None = None,
    ) -> str:
        if self.annotation_dir_name_resolver is None:
            raise DatasetAnnotationPathError(
                "Annotation directory resolver is not configured."
            )

        return self.annotation_dir_name_resolver(dataset_name)

    def get_annotation_file_prefix(
            self,
            *,
            dataset_name: str,
            annotation_dir_name: str,
            group_type: str | None = None,
    ) -> str:
        return annotation_dir_name

    def resolve_annotation_context(self, request) -> dict:
        dataset_name = get_dataset_query_name(request)
        group_by = self.get_group_by(request)
        group_type = self.get_group_type(request)

        annotation_dir_name = self.get_annotation_dir_name(
            dataset_name=dataset_name,
            group_type=group_type,
        )

        annotation_dir = resolve_dataset_annotation_dir(
            annotation_root_dir=self.get_annotation_root_dir(),
            annotation_dir_name=annotation_dir_name,
        )

        file_prefix = self.get_annotation_file_prefix(
            dataset_name=dataset_name,
            annotation_dir_name=annotation_dir_name,
            group_type=group_type,
        )

        survival_file = get_dataset_annotation_survival_file_path(
            annotation_dir=annotation_dir,
            file_prefix=file_prefix,
        )

        return {
            "dataset_name": dataset_name,
            "group_by": group_by,
            "group_type": group_type,
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
            "group_by": context.get("group_by"),
            "group_type": context.get("group_type"),
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

            if (
                self.requires_annotation_dir_name_resolver
                and not self.annotation_dir_name_resolver
            ):
                raise RuntimeError(
                    "Missing annotation_dir_name_resolver."
                )

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

    def get_group_type(self, request):
        return get_timedb_group_type_query(request)

    def get_annotation_dir_name(
        self,
        *,
        dataset_name: str,
        group_type: str | None = None,
    ) -> str:
        return resolve_timedb_group_annotation_dir_name(
            dataset_name=dataset_name,
            group_type=group_type,
        )

    def get_annotation_file_prefix(
        self,
        *,
        dataset_name: str,
        annotation_dir_name: str,
        group_type: str | None = None,
    ) -> str:
        return resolve_timedb_group_annotation_file_prefix(
            dataset_name=dataset_name,
            group_type=group_type,
        )


class SCSTDatasetAnnotationSurvivalKMDataView(
    BaseDatasetAnnotationSurvivalKMDataView
):
    """
    SC/ST Dataset Annotation Kaplan-Meier survival data.

    Query params:
        dataset:
            Dataset name.

        data_type:
            sc | st

        group_by:
            Configured SC/ST Dataset Annotation Group By.

        group_value:
            Selected Group Value.

    Input filename:
        {dataset_name}_survival_analysis_{group_value}.csv

    Directory:
        sc:
            TISCH2_DATASET_ANNOTATIONS_DIR/
                {dataset_name}_{normalized_group_by}/

        st:
            SCTML_DATASET_ANNOTATIONS_DIR/
                {dataset_name}_{normalized_group_by}/

    Source semantics:
        SC/ST Hybrid Reference annotation output.

    KM semantics:
        Reuse the existing common Survival KM builder.

        time:
            n_os

        event:
            c_os_status

        cluster:
            ceRNA_cluster

        groups:
            Cluster_1
            Cluster_2
    """

    source = "SCST"
    network_source_task_type = (
        "SCSTHybridReferenceTask"
    )

    requires_annotation_dir_name_resolver = False

    title = (
        "TCGA-based ceRNA axis survival analysis"
    )

    valid_groups = (
        SCST_HYBRID_REFERENCE_SURVIVAL_GROUPS
    )

    def resolve_annotation_context(
        self,
        request,
    ) -> dict:
        dataset_name = (
            get_dataset_query_name(
                request
            )
        )

        data_type = (
            get_scst_data_type_query(
                request
            )
        )

        group_by = (
            get_scst_group_by_query(
                request
            )
        )

        group_value = (
            get_scst_group_value_query(
                request
            )
        )

        group_by = (
            validate_scst_dataset_group_by(
                dataset_name=dataset_name,
                group_by=group_by,
            )
        )

        annotation_dir = (
            resolve_scst_dataset_group_annotation_dir(
                dataset_name=dataset_name,
                group_by=group_by,
                data_type=data_type,
            )
        )

        survival_file = (
            get_scst_dataset_survival_file_path(
                annotation_dir=annotation_dir,
                dataset_name=dataset_name,
                group_value=group_value,
            )
        )

        return {
            "dataset_name": dataset_name,
            "data_type": data_type,

            "group_by": group_by,
            "group_type": None,
            "group_value": group_value,

            "annotation_dir_name": (
                annotation_dir.name
            ),

            # SC/ST filenames are group-value based.
            # Keep the dataset name as the informational prefix.
            "annotation_file_prefix": (
                dataset_name
            ),

            "annotation_dir": (
                annotation_dir
            ),

            "survival_file": (
                survival_file
            ),
        }

    def get_base_response(
        self,
        context: dict,
    ) -> dict:
        response_data = (
            super().get_base_response(
                context
            )
        )

        response_data.update(
            {
                "data_type": (
                    context.get(
                        "data_type"
                    )
                ),
                "group_value": (
                    context.get(
                        "group_value"
                    )
                ),
            }
        )

        return response_data
