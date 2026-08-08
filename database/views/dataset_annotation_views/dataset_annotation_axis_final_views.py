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
    PAIRED_COHORT_AXIS_FINAL_REQUIRED_COLUMNS,
    HYBRID_REFERENCE_AXIS_FINAL_COLUMNS,
    HYBRID_REFERENCE_AXIS_FINAL_REQUIRED_COLUMNS,
    SCST_HYBRID_REFERENCE_AXIS_FINAL_COLUMNS,
    SCST_HYBRID_REFERENCE_AXIS_FINAL_REQUIRED_COLUMNS,
    WORKFLOW_AXIS_FINAL_NUMERIC_COLUMNS,
)
from database.utils.dataset_annotation_utils.path_utils import (
    DatasetAnnotationInputError,
    DatasetAnnotationPathError,
    get_dataset_query_name,
    resolve_dataset_annotation_dir,
    get_dataset_annotation_axis_final_file_path,
    resolve_tcga_annotation_dir_name,
    resolve_timedb_annotation_dir_name,
    get_timedb_group_type_query,
    resolve_timedb_group_annotation_dir_name,
    resolve_timedb_group_annotation_file_prefix,
)
from database.utils.dataset_annotation_utils.scst_path_utils import (
    get_scst_data_type_query,
    get_scst_group_by_query,
    get_scst_group_value_query,
    validate_scst_dataset_group_by,
    resolve_scst_dataset_group_annotation_dir,
    get_scst_dataset_axis_final_file_path,
)


class BaseDatasetAnnotationAxisFinalDataView(APIView):
    """
    Shared Dataset Annotation Axis Final endpoint.

    Source-specific views only need to resolve the request context
    and the Axis Final input file. CSV parsing and response
    serialization continue to use the workflow Axis Final utilities.
    """

    source = None
    network_source_task_type = None
    annotation_root_setting_name = None
    annotation_dir_name_resolver = None

    axis_final_columns = []
    axis_final_required_columns = []
    axis_final_numeric_columns = (
        WORKFLOW_AXIS_FINAL_NUMERIC_COLUMNS
    )

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
                f"{self.annotation_root_setting_name} "
                "is not configured."
            )

        return annotation_root_dir

    def get_group_type(
        self,
        request,
    ):
        return None

    def get_group_by(
        self,
        request,
    ):
        group_by = request.query_params.get(
            "group_by"
        )

        if group_by is None:
            return None

        group_by = str(
            group_by
        ).strip()

        return group_by or None

    def get_data_type(
        self,
        request,
    ):
        """
        Optional source-specific data type.

        TCGA / TIMEDB:
            None

        SC/ST:
            sc | st
        """
        return None

    def get_group_value(
        self,
        request,
    ):
        """
        Optional visualization-level group value.

        TCGA / TIMEDB:
            None

        SC/ST:
            required by the source-specific view.
        """
        return None

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

        return self.annotation_dir_name_resolver(
            dataset_name
        )

    def get_annotation_file_prefix(
        self,
        *,
        dataset_name: str,
        annotation_dir_name: str,
        group_type: str | None = None,
    ) -> str:
        return annotation_dir_name

    def get_axis_final_file_info(
        self,
        *,
        dataset_name: str,
        group_by: str | None = None,
        group_type: str | None = None,
        data_type: str | None = None,
        group_value: str | None = None,
    ) -> dict:
        """
        Default TCGA/TIMEDB resolver.

        The expanded context arguments allow SC/ST to reuse the same
        Base request pipeline without overriding get().
        """
        annotation_dir_name = (
            self.get_annotation_dir_name(
                dataset_name=dataset_name,
                group_type=group_type,
            )
        )

        annotation_dir = (
            resolve_dataset_annotation_dir(
                annotation_root_dir=(
                    self.get_annotation_root_dir()
                ),
                annotation_dir_name=(
                    annotation_dir_name
                ),
            )
        )

        file_prefix = (
            self.get_annotation_file_prefix(
                dataset_name=dataset_name,
                annotation_dir_name=(
                    annotation_dir_name
                ),
                group_type=group_type,
            )
        )

        axis_final_file = (
            get_dataset_annotation_axis_final_file_path(
                annotation_dir=annotation_dir,
                file_prefix=file_prefix,
            )
        )

        return {
            "dataset_name": dataset_name,
            "group_type": group_type,
            "annotation_dir_name": (
                annotation_dir_name
            ),
            "annotation_file_prefix": (
                file_prefix
            ),
            "axis_final_file": (
                axis_final_file
            ),
        }

    def get(self, request):
        try:
            if not self.source:
                raise RuntimeError(
                    "Missing source."
                )

            if not self.network_source_task_type:
                raise RuntimeError(
                    "Missing network_source_task_type."
                )

            if not self.axis_final_columns:
                raise RuntimeError(
                    "Missing axis_final_columns."
                )

            dataset_name = (
                get_dataset_query_name(
                    request
                )
            )

            group_by = (
                self.get_group_by(
                    request
                )
            )

            group_type = (
                self.get_group_type(
                    request
                )
            )

            data_type = (
                self.get_data_type(
                    request
                )
            )

            group_value = (
                self.get_group_value(
                    request
                )
            )

            required_columns = (
                self.axis_final_required_columns
                or self.axis_final_columns
            )

            file_info = (
                self.get_axis_final_file_info(
                    dataset_name=dataset_name,
                    group_by=group_by,
                    group_type=group_type,
                    data_type=data_type,
                    group_value=group_value,
                )
            )

            axis_file, df = (
                read_axis_final_file_by_path(
                    file_path=(
                        file_info[
                            "axis_final_file"
                        ]
                    ),
                    required_columns=(
                        required_columns
                    ),
                )
            )

            base_response = {
                "success": True,
                "source": self.source,
                "dataset_name": (
                    file_info[
                        "dataset_name"
                    ]
                ),
                "data_type": data_type,
                "group_by": group_by,
                "group_type": group_type,
                "group_value": group_value,
                "annotation_dir_name": (
                    file_info[
                        "annotation_dir_name"
                    ]
                ),
                "annotation_file_prefix": (
                    file_info[
                        "annotation_file_prefix"
                    ]
                ),
                "network_source_task_type": (
                    self.network_source_task_type
                ),
            }

            response_data = (
                build_axis_final_response_from_dataframe(
                    df=df,
                    axis_file_name=(
                        axis_file.name
                    ),
                    columns=(
                        self.axis_final_columns
                    ),
                    required_columns=(
                        required_columns
                    ),
                    numeric_columns=(
                        self.axis_final_numeric_columns
                    ),
                    base_response=(
                        base_response
                    ),
                )
            )

            return Response(
                response_data,
                status=status.HTTP_200_OK,
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
            ValueError,
        ) as e:
            return Response(
                {
                    "success": False,
                    "detail": str(e),
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        except Exception as e:
            print(
                traceback.format_exc()
            )

            return Response(
                {
                    "success": False,
                    "detail": (
                        f"Server error: {str(e)}"
                    ),
                },
                status=(
                    status.HTTP_500_INTERNAL_SERVER_ERROR
                ),
            )


class TCGADatasetAnnotationAxisFinalDataView(
    BaseDatasetAnnotationAxisFinalDataView
):
    """
    TCGA dataset annotation Axis Final data.

    Input:
        ?dataset=TCGA_ACC_mRNA

    Resolution:
        TCGA_ACC_mRNA -> TCGA_ACC

    Source semantics:
        Paired Cohort annotation output.
    """

    source = "TCGA"
    network_source_task_type = (
        "PairedCohortTask"
    )

    annotation_root_setting_name = (
        "TCGA_DATASET_ANNOTATIONS_DIR"
    )
    annotation_dir_name_resolver = staticmethod(
        resolve_tcga_annotation_dir_name
    )

    axis_final_columns = (
        PAIRED_COHORT_AXIS_FINAL_COLUMNS
    )
    axis_final_required_columns = (
        PAIRED_COHORT_AXIS_FINAL_REQUIRED_COLUMNS
    )


class TIMEDBDatasetAnnotationAxisFinalDataView(
    BaseDatasetAnnotationAxisFinalDataView
):
    """
    TIMEDB dataset annotation Axis Final data.

    Source semantics:
        Hybrid Reference annotation output.
    """

    source = "TIMEDB"
    network_source_task_type = (
        "HybridReferenceTask"
    )

    annotation_root_setting_name = (
        "TIMEDB_DATASET_ANNOTATIONS_DIR"
    )
    annotation_dir_name_resolver = staticmethod(
        resolve_timedb_annotation_dir_name
    )

    axis_final_columns = (
        HYBRID_REFERENCE_AXIS_FINAL_COLUMNS
    )
    axis_final_required_columns = (
        HYBRID_REFERENCE_AXIS_FINAL_REQUIRED_COLUMNS
    )

    def get_group_type(
        self,
        request,
    ):
        return (
            get_timedb_group_type_query(
                request
            )
        )

    def get_annotation_dir_name(
        self,
        *,
        dataset_name: str,
        group_type: str | None = None,
    ) -> str:
        return (
            resolve_timedb_group_annotation_dir_name(
                dataset_name=dataset_name,
                group_type=group_type,
            )
        )

    def get_annotation_file_prefix(
        self,
        *,
        dataset_name: str,
        annotation_dir_name: str,
        group_type: str | None = None,
    ) -> str:
        return (
            resolve_timedb_group_annotation_file_prefix(
                dataset_name=dataset_name,
                group_type=group_type,
            )
        )


class SCSTDatasetAnnotationAxisFinalDataView(
    BaseDatasetAnnotationAxisFinalDataView
):
    """
    SC/ST Dataset Annotation Axis Final data.

    Input:
        ?dataset=BCC_GSE123813_aPD1
        &data_type=sc
        &group_by=Celltype major lineage
        &group_value=B

    Directory:
        sc:
            TISCH2_DATASET_ANNOTATIONS_DIR/
                {dataset}_{normalized_group_by}/

        st:
            SCTML_DATASET_ANNOTATIONS_DIR/
                {dataset}_{normalized_group_by}/

    File:
        {dataset}_ceRNA_{group_value}_axis_final.csv

    Source semantics:
        SC/ST Hybrid Reference annotation output.
    """

    source = "SCST"
    network_source_task_type = (
        "SCSTHybridReferenceTask"
    )

    axis_final_columns = (
        SCST_HYBRID_REFERENCE_AXIS_FINAL_COLUMNS
    )
    axis_final_required_columns = (
        SCST_HYBRID_REFERENCE_AXIS_FINAL_REQUIRED_COLUMNS
    )

    def get_data_type(
        self,
        request,
    ):
        return (
            get_scst_data_type_query(
                request
            )
        )

    def get_group_by(
        self,
        request,
    ):
        return (
            get_scst_group_by_query(
                request
            )
        )

    def get_group_value(
        self,
        request,
    ):
        return (
            get_scst_group_value_query(
                request
            )
        )

    def get_axis_final_file_info(
        self,
        *,
        dataset_name: str,
        group_by: str | None = None,
        group_type: str | None = None,
        data_type: str | None = None,
        group_value: str | None = None,
    ) -> dict:
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

        axis_final_file = (
            get_scst_dataset_axis_final_file_path(
                annotation_dir=annotation_dir,
                dataset_name=dataset_name,
                group_value=group_value,
            )
        )

        return {
            "dataset_name": dataset_name,
            "group_type": None,
            "annotation_dir_name": (
                annotation_dir.name
            ),

            # SC/ST filenames insert group_value between
            # dataset name and the Axis Final suffix.
            # Keep this field only as informational metadata.
            "annotation_file_prefix": (
                dataset_name
            ),

            "axis_final_file": (
                axis_final_file
            ),
        }
