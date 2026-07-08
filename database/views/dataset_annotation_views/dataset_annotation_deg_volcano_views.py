import traceback

from django.conf import settings

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.utils.workflow_detail_utils.workflow_deg_volcano_utils import (
    WORKFLOW_DEG_SCOPE_ALL,
    WORKFLOW_DEG_PAIRED_COHORT_RNA_TYPES,
    WORKFLOW_DEG_HYBRID_REFERENCE_RNA_TYPES,
    WORKFLOW_DEG_PAIRED_COHORT_SCOPES,
    WORKFLOW_DEG_HYBRID_REFERENCE_SCOPES,
    WorkflowDEGVolcanoInputError,
    read_deg_file_by_path,
    build_deg_volcano_response_data_from_dataframe,
)
from database.utils.dataset_annotation_utils.path_utils import (
    DatasetAnnotationInputError,
    DatasetAnnotationPathError,
    get_dataset_query_name,
    resolve_dataset_annotation_dir,
    resolve_tcga_annotation_dir_name,
    resolve_timedb_annotation_dir_name,
    get_dataset_annotation_deg_file_path,
    get_dataset_annotation_available_deg_rna_types,
    get_dataset_annotation_available_deg_scopes, get_timedb_group_type_query, resolve_timedb_group_annotation_dir_name,
    resolve_timedb_group_annotation_file_prefix,
)


def parse_bool_query_param(value, default=True) -> bool:
    if value is None:
        return default

    normalized = str(value).strip().lower()

    if normalized in {"1", "true", "yes", "y"}:
        return True

    if normalized in {"0", "false", "no", "n"}:
        return False

    return default


class BaseDatasetAnnotationDEGVolcanoView(APIView):
    """
    Return dataset-level annotation DEG volcano plot data.

    Query params:
        dataset: dataset name
        rna_type: optional
        deg_scope: optional, default = all
        deg_method: optional, default = limma
        use_padj: optional, default = true

    Filename:
        all:
            {annotation_file_prefix}_{deg_method}_{rna_type}.csv

        intersect:
            {annotation_file_prefix}_{deg_method}_{rna_type}_intersect.csv
    """

    source = None
    network_source_task_type = None

    annotation_root_setting_name = None
    annotation_dir_name_resolver = None

    valid_rna_types = []
    valid_deg_scopes = [WORKFLOW_DEG_SCOPE_ALL]

    default_rna_type = None
    default_deg_scope = WORKFLOW_DEG_SCOPE_ALL
    default_deg_method = "limma"
    default_use_padj = False

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

    def resolve_annotation_context(
            self,
            *,
            dataset_name: str,
            group_type: str | None = None,
    ) -> dict:
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

        return {
            "dataset_name": dataset_name,
            "group_type": group_type,
            "annotation_dir_name": annotation_dir_name,
            "annotation_file_prefix": file_prefix,
            "annotation_dir": annotation_dir,
        }

    def get(self, request):
        try:
            if not self.source:
                raise RuntimeError("Missing source.")

            if not self.network_source_task_type:
                raise RuntimeError("Missing network_source_task_type.")

            if not self.valid_rna_types:
                raise RuntimeError("Missing valid_rna_types.")

            dataset_name = get_dataset_query_name(request)
            group_by = self.get_group_by(request)
            group_type = self.get_group_type(request)

            rna_type = str(
                request.query_params.get(
                    "rna_type",
                    self.default_rna_type or "",
                )
            ).strip()

            deg_scope = str(
                request.query_params.get(
                    "deg_scope",
                    self.default_deg_scope,
                )
            ).strip()

            deg_method = str(
                request.query_params.get(
                    "deg_method",
                    self.default_deg_method,
                )
            ).strip()

            use_padj = parse_bool_query_param(
                request.query_params.get("use_padj"),
                default=self.default_use_padj,
            )

            if not rna_type:
                return Response(
                    {
                        "detail": "Missing query parameter: rna_type."
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if rna_type not in self.valid_rna_types:
                return Response(
                    {
                        "detail": (
                            "Invalid rna_type. Allowed values are: "
                            f"{', '.join(self.valid_rna_types)}."
                        ),
                        "valid_rna_types": self.valid_rna_types,
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if deg_scope not in self.valid_deg_scopes:
                return Response(
                    {
                        "detail": (
                            "Invalid deg_scope. Allowed values are: "
                            f"{', '.join(self.valid_deg_scopes)}."
                        ),
                        "valid_deg_scopes": self.valid_deg_scopes,
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not deg_method:
                return Response(
                    {
                        "detail": "Missing query parameter: deg_method."
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                context = self.resolve_annotation_context(
                    dataset_name=dataset_name,
                    group_type=group_type,
                )

                available_deg_rna_types = (
                    get_dataset_annotation_available_deg_rna_types(
                        annotation_dir=context["annotation_dir"],
                        file_prefix=context["annotation_file_prefix"],
                        deg_method=deg_method,
                        valid_rna_types=self.valid_rna_types,
                        deg_scope=WORKFLOW_DEG_SCOPE_ALL,
                    )
                )

                available_deg_scopes = (
                    get_dataset_annotation_available_deg_scopes(
                        annotation_dir=context["annotation_dir"],
                        file_prefix=context["annotation_file_prefix"],
                        deg_method=deg_method,
                        rna_type=rna_type,
                        valid_scopes=self.valid_deg_scopes,
                    )
                )

                deg_file_path = get_dataset_annotation_deg_file_path(
                    annotation_dir=context["annotation_dir"],
                    file_prefix=context["annotation_file_prefix"],
                    deg_method=deg_method,
                    rna_type=rna_type,
                    deg_scope=deg_scope,
                )

                deg_file, df = read_deg_file_by_path(deg_file_path)

            except FileNotFoundError as e:
                return Response(
                    {
                        "detail": str(e),
                        "available_deg_rna_types": (
                            available_deg_rna_types
                            if "available_deg_rna_types" in locals()
                            else []
                        ),
                        "available_deg_scopes": (
                            available_deg_scopes
                            if "available_deg_scopes" in locals()
                            else []
                        ),
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            except (
                DatasetAnnotationInputError,
                DatasetAnnotationPathError,
                WorkflowDEGVolcanoInputError,
            ) as e:
                return Response(
                    {
                        "detail": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            base_response = {
                "success": True,
                "source": self.source,
                "dataset_name": context["dataset_name"],
                "group_by": group_by,
                "group_type": group_type,
                "annotation_dir_name": context["annotation_dir_name"],
                "annotation_file_prefix": context[
                    "annotation_file_prefix"
                ],
                "network_source_task_type": self.network_source_task_type,
            }

            try:
                response_data = build_deg_volcano_response_data_from_dataframe(
                    df=df,
                    deg_file_name=deg_file.name,
                    rna_type=rna_type,
                    deg_scope=deg_scope,
                    deg_method=deg_method,
                    use_padj=use_padj,
                    base_response=base_response,
                )
            except WorkflowDEGVolcanoInputError as e:
                return Response(
                    {
                        "detail": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            response_data.update(
                {
                    "available_deg_rna_types": available_deg_rna_types,
                    "available_deg_scopes": available_deg_scopes,
                }
            )

            return Response(
                response_data,
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


class TCGADatasetAnnotationDEGVolcanoView(
    BaseDatasetAnnotationDEGVolcanoView
):
    """
    TCGA dataset annotation DEG volcano.

    Source semantics:
        Paired Cohort annotation output.
    """

    source = "TCGA"
    network_source_task_type = "PairedCohortTask"

    annotation_root_setting_name = "TCGA_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_tcga_annotation_dir_name
    )

    valid_rna_types = WORKFLOW_DEG_PAIRED_COHORT_RNA_TYPES
    valid_deg_scopes = WORKFLOW_DEG_PAIRED_COHORT_SCOPES

    default_rna_type = "mRNA"
    default_deg_scope = WORKFLOW_DEG_SCOPE_ALL
    default_deg_method = "limma"
    default_use_padj = False


class TIMEDBDatasetAnnotationDEGVolcanoView(
    BaseDatasetAnnotationDEGVolcanoView
):
    """
    TIMEDB dataset annotation DEG volcano.

    Source semantics:
        Hybrid Reference annotation output.
    """

    source = "TIMEDB"
    network_source_task_type = "HybridReferenceTask"

    annotation_root_setting_name = "TIMEDB_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_timedb_annotation_dir_name
    )

    valid_rna_types = WORKFLOW_DEG_HYBRID_REFERENCE_RNA_TYPES
    valid_deg_scopes = WORKFLOW_DEG_HYBRID_REFERENCE_SCOPES

    default_rna_type = "mRNA"
    default_deg_scope = WORKFLOW_DEG_SCOPE_ALL
    default_deg_method = "limma"
    default_use_padj = False

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
