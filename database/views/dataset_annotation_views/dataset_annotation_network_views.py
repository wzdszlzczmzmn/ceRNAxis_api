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
    get_dataset_annotation_file_path,
    get_timedb_group_type_query,
    resolve_dataset_annotation_dir,
    resolve_tcga_annotation_dir_name,
    resolve_timedb_annotation_dir_name,
    resolve_timedb_group_annotation_dir_name,
    resolve_timedb_group_annotation_file_prefix,
)
from database.utils.dataset_annotation_utils.scst_path_utils import (
    get_scst_data_type_query,
    get_scst_group_by_query,
    get_scst_group_value_query,
    validate_scst_dataset_group_by,
    resolve_scst_dataset_group_annotation_dir,
    get_scst_dataset_cerna_axis_file_path,
    get_scst_dataset_map_immune_axis_file_path,
)


DATASET_ANNOTATION_CERNA_AXIS_SUFFIX = "_ceRNA_axis.csv"
DATASET_ANNOTATION_IMMUNE_AXIS_SUFFIX = "_map_immune_axis.csv"


class BaseDatasetAnnotationNetworkView(
    WorkflowOutputNetworkBuilderMixin,
    APIView,
):
    """
    Shared Dataset Annotation network view.

    Source-specific views are responsible only for resolving the
    request context and the two network input files.

    The actual CSV parsing and graph construction are delegated to
    WorkflowOutputNetworkBuilderMixin.
    """

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
        Default rule:
            annotation file prefix == annotation directory name.

        Example:
            TCGA_ACC/
                TCGA_ACC_ceRNA_axis.csv
                TCGA_ACC_map_immune_axis.csv
        """
        return annotation_dir_name

    def get_group_type(self, request):
        return None

    def get_group_by(self, request):
        group_by = request.query_params.get(
            "group_by"
        )

        if group_by is None:
            return None

        group_by = str(
            group_by
        ).strip()

        return group_by or None

    def get_data_type(self, request):
        """
        Optional source-specific data type.

        TCGA / TIMEDB:
            None

        SC/ST:
            sc | st
        """
        return None

    def get_group_value(self, request):
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

    def get_network_files(
        self,
        *,
        dataset_name: str,
        group_by: str | None = None,
        group_type: str | None = None,
        data_type: str | None = None,
        group_value: str | None = None,
    ) -> dict:
        """
        Default TCGA/TIMEDB file resolver.

        The additional group_by / data_type / group_value arguments
        are accepted so source-specific views can share one Base
        request pipeline without overriding get().
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
            )
        )

        ceRNA_axis_file = (
            get_dataset_annotation_file_path(
                annotation_dir=annotation_dir,
                file_prefix=file_prefix,
                filename_suffix=(
                    DATASET_ANNOTATION_CERNA_AXIS_SUFFIX
                ),
            )
        )

        immune_axis_file = (
            get_dataset_annotation_file_path(
                annotation_dir=annotation_dir,
                file_prefix=file_prefix,
                filename_suffix=(
                    DATASET_ANNOTATION_IMMUNE_AXIS_SUFFIX
                ),
            )
        )

        return {
            "dataset_name": dataset_name,
            "annotation_dir_name": (
                annotation_dir_name
            ),
            "annotation_file_prefix": (
                file_prefix
            ),
            "annotation_dir": annotation_dir,
            "ceRNA_axis_file": (
                ceRNA_axis_file
            ),
            "immune_axis_file": (
                immune_axis_file
            ),
        }

    def get(self, request, *args, **kwargs):
        try:
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

            file_info = (
                self.get_network_files(
                    dataset_name=dataset_name,
                    group_by=group_by,
                    group_type=group_type,
                    data_type=data_type,
                    group_value=group_value,
                )
            )

            result = (
                self.build_network_from_files(
                    ceRNA_axis_file=(
                        file_info[
                            "ceRNA_axis_file"
                        ]
                    ),
                    immune_axis_file=(
                        file_info[
                            "immune_axis_file"
                        ]
                    ),
                )
            )

            result.update(
                {
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
            )

            return Response(
                result,
                status=status.HTTP_200_OK,
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

        except Exception as e:
            print(
                traceback.format_exc()
            )

            return Response(
                {
                    "success": False,
                    "msg": (
                        f"Server error: {str(e)}"
                    ),
                },
                status=(
                    status.HTTP_500_INTERNAL_SERVER_ERROR
                ),
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
    network_source_task_type = (
        "PairedCohortTask"
    )
    annotation_root_setting_name = (
        "TCGA_DATASET_ANNOTATIONS_DIR"
    )
    annotation_dir_name_resolver = staticmethod(
        resolve_tcga_annotation_dir_name
    )

    cerna_node_source = (
        "paired_cohort_cerna_axis"
    )
    cerna_edge_source = (
        "paired_cohort_cerna_axis"
    )
    immune_node_source = (
        "paired_cohort_immune_axis"
    )
    immune_edge_source = (
        "paired_cohort_immune_axis"
    )


class TIMEDBDatasetAnnotationNetworkView(
    BaseDatasetAnnotationNetworkView
):
    """
    TIMEDB dataset annotation network.

    Input:
        ?dataset=GSE19750
        &group_by=grade
        &group_type=grade

    Annotation source:
        Hybrid Reference workflow output.
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

    cerna_node_source = (
        "hybrid_reference_cerna_axis"
    )
    cerna_edge_source = (
        "hybrid_reference_cerna_axis"
    )
    immune_node_source = (
        "hybrid_reference_immune_axis"
    )
    immune_edge_source = (
        "hybrid_reference_immune_axis"
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
    ) -> str:
        """
        TIMEDB grouped annotations use different directories,
        while the file prefix remains the original dataset name.

        Examples:
            GSE20194/
                GSE20194_ceRNA_axis.csv

            GSE20194_grade/
                GSE20194_ceRNA_axis.csv

            GSE20194_stage/
                GSE20194_ceRNA_axis.csv
        """
        return (
            resolve_timedb_group_annotation_file_prefix(
                dataset_name=dataset_name,
                group_type=None,
            )
        )


class SCSTDatasetAnnotationNetworkView(
    BaseDatasetAnnotationNetworkView
):
    """
    SC/ST dataset annotation network.

    Input:
        ?dataset=BCC_GSE123813_aPD1
        &data_type=sc
        &group_by=Celltype major lineage
        &group_value=B

    Directory resolution:
        sc:
            TISCH2_DATASET_ANNOTATIONS_DIR/
                {dataset}_{normalized_group_by}/

        st:
            SCTML_DATASET_ANNOTATIONS_DIR/
                {dataset}_{normalized_group_by}/

    File expectation inside the selected group-by directory:
        {dataset}_ceRNA_{group_value}_axis.csv
        {dataset}_map_immune_axis_{group_value}.csv

    Annotation source:
        SC/ST Hybrid Reference workflow output.
    """

    source = "SCST"
    network_source_task_type = (
        "SCSTHybridReferenceTask"
    )

    cerna_node_source = (
        "scst_hybrid_reference_cerna_axis"
    )
    cerna_edge_source = (
        "scst_hybrid_reference_cerna_axis"
    )
    immune_node_source = (
        "scst_hybrid_reference_immune_axis"
    )
    immune_edge_source = (
        "scst_hybrid_reference_immune_axis"
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

    def get_network_files(
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

        ceRNA_axis_file = (
            get_scst_dataset_cerna_axis_file_path(
                annotation_dir=annotation_dir,
                dataset_name=dataset_name,
                group_value=group_value,
            )
        )

        immune_axis_file = (
            get_scst_dataset_map_immune_axis_file_path(
                annotation_dir=annotation_dir,
                dataset_name=dataset_name,
                group_value=group_value,
            )
        )

        return {
            "dataset_name": dataset_name,
            "annotation_dir_name": (
                annotation_dir.name
            ),
            # SC/ST filenames contain group_value in the middle,
            # so annotation_file_prefix is informational only.
            "annotation_file_prefix": (
                dataset_name
            ),
            "annotation_dir": annotation_dir,
            "ceRNA_axis_file": (
                ceRNA_axis_file
            ),
            "immune_axis_file": (
                immune_axis_file
            ),
        }
