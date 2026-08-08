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
    SCST_WORKFLOW_DEG_RNA_TYPES,
    SCST_WORKFLOW_DEG_SCOPES,
    SCST_WORKFLOW_DEG_SCOPE_ALL,
    SCST_WORKFLOW_DEG_SCOPE_INTERSECT,
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
    get_dataset_annotation_available_deg_scopes,
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
    get_scst_dataset_deg_file_path,
    get_scst_dataset_deg_intersect_file_path,
    is_scst_existing_file,
)


def parse_bool_query_param(
    value,
    default=True,
) -> bool:
    if value is None:
        return default

    normalized = str(
        value
    ).strip().lower()

    if normalized in {
        "1",
        "true",
        "yes",
        "y",
    }:
        return True

    if normalized in {
        "0",
        "false",
        "no",
        "n",
    }:
        return False

    return default


class BaseDatasetAnnotationDEGVolcanoView(APIView):
    """
    Shared Dataset Annotation DEG Volcano endpoint.

    Query params:
        dataset
        rna_type
        deg_scope
        deg_method
        use_padj

    Source-specific implementations may additionally require:
        data_type
        group_by
        group_value

    TCGA/TIMEDB use the generic Dataset Annotation DEG filename
    helpers. SC/ST overrides the DEG file/availability hooks because
    its files are group-value based:

        all:
            {dataset}_deg_{group_value}.csv

        intersect:
            {dataset}_mRNA_deg_{group_value}_intersect.csv
    """

    source = None
    network_source_task_type = None

    annotation_root_setting_name = None
    annotation_dir_name_resolver = None

    valid_rna_types = []
    valid_deg_scopes = [
        WORKFLOW_DEG_SCOPE_ALL,
    ]

    default_rna_type = None
    default_deg_scope = (
        WORKFLOW_DEG_SCOPE_ALL
    )
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
        TCGA / TIMEDB:
            None

        SC/ST:
            required.
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

    def resolve_annotation_context(
        self,
        *,
        dataset_name: str,
        group_by: str | None = None,
        group_type: str | None = None,
        data_type: str | None = None,
        group_value: str | None = None,
    ) -> dict:
        """
        Default TCGA/TIMEDB annotation context.
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

        return {
            "dataset_name": dataset_name,
            "group_type": group_type,
            "annotation_dir_name": (
                annotation_dir_name
            ),
            "annotation_file_prefix": (
                file_prefix
            ),
            "annotation_dir": (
                annotation_dir
            ),
        }

    def get_available_deg_rna_types(
        self,
        *,
        context: dict,
        deg_method: str,
        group_value: str | None,
    ) -> list[str]:
        """
        Default TCGA/TIMEDB DEG RNA-type availability.
        """
        return (
            get_dataset_annotation_available_deg_rna_types(
                annotation_dir=(
                    context[
                        "annotation_dir"
                    ]
                ),
                file_prefix=(
                    context[
                        "annotation_file_prefix"
                    ]
                ),
                deg_method=deg_method,
                valid_rna_types=(
                    self.valid_rna_types
                ),
                deg_scope=(
                    WORKFLOW_DEG_SCOPE_ALL
                ),
            )
        )

    def get_available_deg_scopes(
        self,
        *,
        context: dict,
        deg_method: str,
        rna_type: str,
        group_value: str | None,
    ) -> list[str]:
        """
        Default TCGA/TIMEDB DEG-scope availability.
        """
        return (
            get_dataset_annotation_available_deg_scopes(
                annotation_dir=(
                    context[
                        "annotation_dir"
                    ]
                ),
                file_prefix=(
                    context[
                        "annotation_file_prefix"
                    ]
                ),
                deg_method=deg_method,
                rna_type=rna_type,
                valid_scopes=(
                    self.valid_deg_scopes
                ),
            )
        )

    def get_deg_file_path(
        self,
        *,
        context: dict,
        deg_method: str,
        rna_type: str,
        deg_scope: str,
        group_value: str | None,
    ):
        """
        Default TCGA/TIMEDB DEG file resolver.
        """
        return (
            get_dataset_annotation_deg_file_path(
                annotation_dir=(
                    context[
                        "annotation_dir"
                    ]
                ),
                file_prefix=(
                    context[
                        "annotation_file_prefix"
                    ]
                ),
                deg_method=deg_method,
                rna_type=rna_type,
                deg_scope=deg_scope,
            )
        )

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

            if not self.valid_rna_types:
                raise RuntimeError(
                    "Missing valid_rna_types."
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
                request.query_params.get(
                    "use_padj"
                ),
                default=(
                    self.default_use_padj
                ),
            )

            if not rna_type:
                return Response(
                    {
                        "detail": (
                            "Missing query parameter: rna_type."
                        ),
                    },
                    status=(
                        status.HTTP_400_BAD_REQUEST
                    ),
                )

            if rna_type not in self.valid_rna_types:
                return Response(
                    {
                        "detail": (
                            "Invalid rna_type. "
                            "Allowed values are: "
                            f"{', '.join(self.valid_rna_types)}."
                        ),
                        "valid_rna_types": (
                            self.valid_rna_types
                        ),
                    },
                    status=(
                        status.HTTP_400_BAD_REQUEST
                    ),
                )

            if deg_scope not in self.valid_deg_scopes:
                return Response(
                    {
                        "detail": (
                            "Invalid deg_scope. "
                            "Allowed values are: "
                            f"{', '.join(self.valid_deg_scopes)}."
                        ),
                        "valid_deg_scopes": (
                            self.valid_deg_scopes
                        ),
                    },
                    status=(
                        status.HTTP_400_BAD_REQUEST
                    ),
                )

            if not deg_method:
                return Response(
                    {
                        "detail": (
                            "Missing query parameter: deg_method."
                        ),
                    },
                    status=(
                        status.HTTP_400_BAD_REQUEST
                    ),
                )

            context = (
                self.resolve_annotation_context(
                    dataset_name=dataset_name,
                    group_by=group_by,
                    group_type=group_type,
                    data_type=data_type,
                    group_value=group_value,
                )
            )

            available_deg_rna_types = (
                self.get_available_deg_rna_types(
                    context=context,
                    deg_method=deg_method,
                    group_value=group_value,
                )
            )

            available_deg_scopes = (
                self.get_available_deg_scopes(
                    context=context,
                    deg_method=deg_method,
                    rna_type=rna_type,
                    group_value=group_value,
                )
            )

            deg_file_path = (
                self.get_deg_file_path(
                    context=context,
                    deg_method=deg_method,
                    rna_type=rna_type,
                    deg_scope=deg_scope,
                    group_value=group_value,
                )
            )

            deg_file, df = (
                read_deg_file_by_path(
                    deg_file_path
                )
            )

            base_response = {
                "success": True,
                "source": self.source,
                "dataset_name": (
                    context[
                        "dataset_name"
                    ]
                ),
                "data_type": data_type,
                "group_by": group_by,
                "group_type": group_type,
                "group_value": group_value,
                "annotation_dir_name": (
                    context[
                        "annotation_dir_name"
                    ]
                ),
                "annotation_file_prefix": (
                    context[
                        "annotation_file_prefix"
                    ]
                ),
                "network_source_task_type": (
                    self.network_source_task_type
                ),
            }

            response_data = (
                build_deg_volcano_response_data_from_dataframe(
                    df=df,
                    deg_file_name=(
                        deg_file.name
                    ),
                    rna_type=rna_type,
                    deg_scope=deg_scope,
                    deg_method=deg_method,
                    use_padj=use_padj,
                    base_response=(
                        base_response
                    ),
                )
            )

            response_data.update(
                {
                    "available_deg_rna_types": (
                        available_deg_rna_types
                    ),
                    "available_deg_scopes": (
                        available_deg_scopes
                    ),
                }
            )

            return Response(
                response_data,
                status=status.HTTP_200_OK,
            )

        except FileNotFoundError as e:
            return Response(
                {
                    "detail": str(e),
                    "available_deg_rna_types": (
                        available_deg_rna_types
                        if "available_deg_rna_types"
                        in locals()
                        else []
                    ),
                    "available_deg_scopes": (
                        available_deg_scopes
                        if "available_deg_scopes"
                        in locals()
                        else []
                    ),
                },
                status=status.HTTP_404_NOT_FOUND,
            )

        except (
            DatasetAnnotationInputError,
            DatasetAnnotationPathError,
            WorkflowDEGVolcanoInputError,
            ValueError,
        ) as e:
            return Response(
                {
                    "detail": str(e),
                },
                status=(
                    status.HTTP_400_BAD_REQUEST
                ),
            )

        except Exception as e:
            print(
                traceback.format_exc()
            )

            return Response(
                {
                    "detail": (
                        f"Server error: {str(e)}"
                    ),
                },
                status=(
                    status.HTTP_500_INTERNAL_SERVER_ERROR
                ),
            )


class TCGADatasetAnnotationDEGVolcanoView(
    BaseDatasetAnnotationDEGVolcanoView
):
    """
    TCGA Dataset Annotation DEG Volcano.

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

    valid_rna_types = (
        WORKFLOW_DEG_PAIRED_COHORT_RNA_TYPES
    )
    valid_deg_scopes = (
        WORKFLOW_DEG_PAIRED_COHORT_SCOPES
    )

    default_rna_type = "mRNA"
    default_deg_scope = (
        WORKFLOW_DEG_SCOPE_ALL
    )
    default_deg_method = "limma"
    default_use_padj = False


class TIMEDBDatasetAnnotationDEGVolcanoView(
    BaseDatasetAnnotationDEGVolcanoView
):
    """
    TIMEDB Dataset Annotation DEG Volcano.

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

    valid_rna_types = (
        WORKFLOW_DEG_HYBRID_REFERENCE_RNA_TYPES
    )
    valid_deg_scopes = (
        WORKFLOW_DEG_HYBRID_REFERENCE_SCOPES
    )

    default_rna_type = "mRNA"
    default_deg_scope = (
        WORKFLOW_DEG_SCOPE_ALL
    )
    default_deg_method = "limma"
    default_use_padj = False

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


class SCSTDatasetAnnotationDEGVolcanoView(
    BaseDatasetAnnotationDEGVolcanoView
):
    """
    SC/ST Dataset Annotation DEG Volcano.

    Input:
        ?dataset=BCC_GSE123813_aPD1
        &data_type=sc
        &group_by=Celltype major lineage
        &group_value=B
        &rna_type=mRNA
        &deg_scope=all
        &deg_method=limma
        &use_padj=false

    Current SC/ST DEG constraints:
        RNA type:
            mRNA

        scopes:
            all
            intersect

    Files:
        all:
            {dataset}_deg_{group_value}.csv

        intersect:
            {dataset}_mRNA_deg_{group_value}_intersect.csv

    Availability rule:
        - all must exist for mRNA DEG to be available.
        - intersect is exposed only when all also exists.
        - an orphan intersect file is not treated as available.
    """

    source = "SCST"
    network_source_task_type = (
        "SCSTHybridReferenceTask"
    )

    valid_rna_types = (
        SCST_WORKFLOW_DEG_RNA_TYPES
    )
    valid_deg_scopes = (
        SCST_WORKFLOW_DEG_SCOPES
    )

    default_rna_type = "mRNA"
    default_deg_scope = (
        SCST_WORKFLOW_DEG_SCOPE_ALL
    )

    # SC/ST filenames do not contain a DEG-method token.
    # Keep limma as Dataset Annotation response/query metadata.
    default_deg_method = "limma"

    # Dataset Annotation uses raw p-value by default.
    default_use_padj = False

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

    def resolve_annotation_context(
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

        return {
            "dataset_name": dataset_name,
            "group_type": None,
            "annotation_dir_name": (
                annotation_dir.name
            ),

            # SC/ST DEG filenames are group-value based and do
            # not use the generic file-prefix convention.
            "annotation_file_prefix": (
                dataset_name
            ),

            "annotation_dir": (
                annotation_dir
            ),
        }

    def get_available_deg_rna_types(
        self,
        *,
        context: dict,
        deg_method: str,
        group_value: str | None,
    ) -> list[str]:
        """
        SC/ST availability is anchored to the all-scope file.

        Truth table:
            all=yes, intersect=no
                -> ["mRNA"]

            all=yes, intersect=yes
                -> ["mRNA"]

            all=no, intersect=yes
                -> []

            all=no, intersect=no
                -> []
        """
        all_file = (
            get_scst_dataset_deg_file_path(
                annotation_dir=(
                    context[
                        "annotation_dir"
                    ]
                ),
                dataset_name=(
                    context[
                        "dataset_name"
                    ]
                ),
                group_value=group_value,
            )
        )

        if not is_scst_existing_file(
            all_file
        ):
            return []

        return [
            "mRNA",
        ]

    def get_available_deg_scopes(
        self,
        *,
        context: dict,
        deg_method: str,
        rna_type: str,
        group_value: str | None,
    ) -> list[str]:
        """
        Do not expose an orphan intersect file.
        """
        all_file = (
            get_scst_dataset_deg_file_path(
                annotation_dir=(
                    context[
                        "annotation_dir"
                    ]
                ),
                dataset_name=(
                    context[
                        "dataset_name"
                    ]
                ),
                group_value=group_value,
            )
        )

        if not is_scst_existing_file(
            all_file
        ):
            return []

        available_scopes = [
            SCST_WORKFLOW_DEG_SCOPE_ALL,
        ]

        intersect_file = (
            get_scst_dataset_deg_intersect_file_path(
                annotation_dir=(
                    context[
                        "annotation_dir"
                    ]
                ),
                dataset_name=(
                    context[
                        "dataset_name"
                    ]
                ),
                group_value=group_value,
            )
        )

        if (
            SCST_WORKFLOW_DEG_SCOPE_INTERSECT
            in self.valid_deg_scopes
            and is_scst_existing_file(
                intersect_file
            )
        ):
            available_scopes.append(
                SCST_WORKFLOW_DEG_SCOPE_INTERSECT
            )

        return available_scopes

    def get_deg_file_path(
        self,
        *,
        context: dict,
        deg_method: str,
        rna_type: str,
        deg_scope: str,
        group_value: str | None,
    ):
        if rna_type != "mRNA":
            raise DatasetAnnotationInputError(
                "SC/ST Dataset Annotation DEG "
                "currently supports mRNA only."
            )

        all_file = (
            get_scst_dataset_deg_file_path(
                annotation_dir=(
                    context[
                        "annotation_dir"
                    ]
                ),
                dataset_name=(
                    context[
                        "dataset_name"
                    ]
                ),
                group_value=group_value,
            )
        )

        if deg_scope == (
            SCST_WORKFLOW_DEG_SCOPE_ALL
        ):
            return all_file

        if deg_scope == (
            SCST_WORKFLOW_DEG_SCOPE_INTERSECT
        ):
            # Enforce the same truth table as the availability API.
            if not is_scst_existing_file(
                all_file
            ):
                raise FileNotFoundError(
                    "SC/ST all-scope DEG file is not "
                    "available, so intersect DEG is not "
                    "considered available."
                )

            return (
                get_scst_dataset_deg_intersect_file_path(
                    annotation_dir=(
                        context[
                            "annotation_dir"
                        ]
                    ),
                    dataset_name=(
                        context[
                            "dataset_name"
                        ]
                    ),
                    group_value=group_value,
                )
            )

        raise DatasetAnnotationInputError(
            f"Invalid SC/ST DEG scope: {deg_scope}."
        )
