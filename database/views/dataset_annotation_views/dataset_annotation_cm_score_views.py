import traceback
from pathlib import Path

from django.conf import settings

from rest_framework import status
from rest_framework.response import Response
from rest_framework.views import APIView

from analysis.utils.workflow_detail_utils.workflow_cm_score_utils import (
    WorkflowCMScoreInputError,
    WorkflowCMScorePathError,
    validate_cm_score_item_value,
)
from database.utils.dataset_annotation_utils.dataset_annotation_cm_score_utils import (
    build_dataset_annotation_cm_score_options_response,
    read_dataset_annotation_cm_score_file,
    build_dataset_annotation_cm_score_result_response,
)
from database.utils.dataset_annotation_utils.path_utils import (
    DatasetAnnotationInputError,
    DatasetAnnotationPathError,
    get_dataset_annotation_cmdrug_result_dir_path,
    get_dataset_query_name,
    get_timedb_group_type_query,
    resolve_dataset_annotation_dir,
    resolve_tcga_annotation_dir_name,
    resolve_timedb_annotation_dir_name,
    resolve_timedb_group_annotation_dir_name,
    resolve_timedb_group_annotation_file_prefix,
    validate_timedb_group_selection,
)
from database.utils.dataset_annotation_utils.scst_path_utils import (
    get_scst_data_type_query,
    get_scst_group_by_query,
    get_scst_group_value_query,
    get_scst_dataset_cmdrug_result_dir_path,
    resolve_scst_dataset_group_annotation_dir,
    validate_scst_dataset_group_by,
)


class BaseDatasetAnnotationCMScoreContextMixin:
    """
    Resolve Dataset Annotation CM-score context.

    TCGA and TIMEDB use the generic annotation directory + file
    prefix convention:

        {annotation_dir}/
            {file_prefix}_CMdrug_result/
                {item}_CM_scores.csv

    SC/ST overrides resolve_annotation_context() because its CMdrug
    directory is group-value-specific:

        {dataset}_{group_by}/
            {dataset}_CMdrug_result_{group_value}/
                {item}_CM_scores.csv
    """

    source = None
    network_source_task_type = None

    annotation_root_setting_name = None
    annotation_dir_name_resolver = None

    uses_generic_annotation_context = True

    def validate_configuration(self) -> None:
        if not self.source:
            raise RuntimeError(
                "Missing Dataset Annotation CM Score source."
            )

        if not self.network_source_task_type:
            raise RuntimeError(
                "Missing network_source_task_type."
            )

        if not self.uses_generic_annotation_context:
            return

        if not self.annotation_root_setting_name:
            raise RuntimeError(
                "Missing annotation_root_setting_name."
            )

        if self.annotation_dir_name_resolver is None:
            raise RuntimeError(
                "Missing annotation_dir_name_resolver."
            )

    def get_annotation_root_dir(self) -> Path:
        if not self.annotation_root_setting_name:
            raise DatasetAnnotationPathError(
                "Annotation root setting name is not configured."
            )

        root_value = getattr(
            settings,
            self.annotation_root_setting_name,
            None,
        )

        if not root_value:
            raise DatasetAnnotationPathError(
                f"{self.annotation_root_setting_name} "
                "is not configured."
            )

        root_dir = Path(
            root_value
        ).resolve()

        if (
            not root_dir.exists()
            or not root_dir.is_dir()
        ):
            raise DatasetAnnotationPathError(
                "Dataset annotation root directory "
                "is not available."
            )

        return root_dir

    def get_group_by(
        self,
        request,
    ) -> str | None:
        group_by = request.query_params.get(
            "group_by"
        )

        if group_by is None:
            return None

        group_by = str(
            group_by
        ).strip()

        return group_by or None

    def get_group_type(
        self,
        request,
    ) -> str | None:
        return None

    def validate_group_context(
        self,
        *,
        dataset_name: str,
        group_by: str | None,
        group_type: str | None,
    ) -> tuple[str | None, str | None]:
        return (
            group_by,
            group_type,
        )

    def get_annotation_dir_name(
        self,
        *,
        dataset_name: str,
        group_type: str | None = None,
    ) -> str:
        if self.annotation_dir_name_resolver is None:
            raise DatasetAnnotationPathError(
                "Annotation directory resolver "
                "is not configured."
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
        request,
    ) -> dict:
        dataset_name = get_dataset_query_name(
            request
        )

        group_by = self.get_group_by(
            request
        )

        group_type = self.get_group_type(
            request
        )

        (
            group_by,
            group_type,
        ) = self.validate_group_context(
            dataset_name=dataset_name,
            group_by=group_by,
            group_type=group_type,
        )

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
                annotation_dir_name=annotation_dir_name,
                group_type=group_type,
            )
        )

        cm_results_dir = (
            get_dataset_annotation_cmdrug_result_dir_path(
                annotation_dir=annotation_dir,
                file_prefix=file_prefix,
            )
        )

        return {
            "dataset_name": dataset_name,

            "group_by": group_by,
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

            "cm_results_dir": (
                cm_results_dir
            ),
        }

    def get_base_response(
        self,
        context: dict,
    ) -> dict:
        response_data = {
            "success": True,

            "source": self.source,

            "dataset_name": (
                context[
                    "dataset_name"
                ]
            ),

            "group_by": (
                context.get(
                    "group_by"
                )
            ),

            "group_type": (
                context.get(
                    "group_type"
                )
            ),

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

        if (
            context.get(
                "data_type"
            )
            is not None
        ):
            response_data[
                "data_type"
            ] = context[
                "data_type"
            ]

        if (
            context.get(
                "group_value"
            )
            is not None
        ):
            response_data[
                "group_value"
            ] = context[
                "group_value"
            ]

        return response_data


class BaseDatasetAnnotationCMScoreView(
    BaseDatasetAnnotationCMScoreContextMixin,
    APIView,
):
    """
    Shared HTTP/error layer for Dataset Annotation CM Score.
    """

    def build_response_data(
        self,
        *,
        request,
        context: dict,
    ) -> dict:
        raise NotImplementedError

    def get(
        self,
        request,
        *args,
        **kwargs,
    ):
        try:
            self.validate_configuration()

            context = (
                self.resolve_annotation_context(
                    request
                )
            )

            response_data = (
                self.build_response_data(
                    request=request,
                    context=context,
                )
            )

            return Response(
                response_data,
                status=status.HTTP_200_OK,
            )

        except (
            DatasetAnnotationInputError,
            WorkflowCMScoreInputError,
        ) as exc:
            return Response(
                {
                    "success": False,
                    "detail": str(
                        exc
                    ),
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        except FileNotFoundError as exc:
            return Response(
                {
                    "success": False,
                    "detail": str(
                        exc
                    ),
                },
                status=status.HTTP_404_NOT_FOUND,
            )

        except (
            DatasetAnnotationPathError,
            WorkflowCMScorePathError,
        ) as exc:
            return Response(
                {
                    "success": False,
                    "detail": str(
                        exc
                    ),
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        except Exception as exc:
            print(
                traceback.format_exc()
            )

            return Response(
                {
                    "success": False,
                    "detail": (
                        f"Server error: {str(exc)}"
                    ),
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )


class BaseDatasetAnnotationCMScoreOptionsView(
    BaseDatasetAnnotationCMScoreView
):
    """
    Return available {item}_CM_scores.csv files.

    A missing CM-results directory is represented as an empty
    result set rather than a 404, matching Workflow behavior.
    """

    def build_response_data(
        self,
        *,
        request,
        context: dict,
    ) -> dict:
        return (
            build_dataset_annotation_cm_score_options_response(
                base_response=(
                    self.get_base_response(
                        context
                    )
                ),
                cm_results_dir=(
                    context[
                        "cm_results_dir"
                    ]
                ),
            )
        )


class BaseDatasetAnnotationCMScoreResultView(
    BaseDatasetAnnotationCMScoreView
):
    """
    Return all rows from one selected {item}_CM_scores.csv file.
    """

    def build_response_data(
        self,
        *,
        request,
        context: dict,
    ) -> dict:
        item_value = (
            validate_cm_score_item_value(
                request.query_params.get(
                    "item"
                )
            )
        )

        (
            file_path,
            dataframe,
        ) = (
            read_dataset_annotation_cm_score_file(
                cm_results_dir=(
                    context[
                        "cm_results_dir"
                    ]
                ),
                item_value=item_value,
            )
        )

        return (
            build_dataset_annotation_cm_score_result_response(
                base_response=(
                    self.get_base_response(
                        context
                    )
                ),
                item_value=item_value,
                file_path=file_path,
                dataframe=dataframe,
            )
        )


class TCGADatasetAnnotationCMScoreMixin:
    """
    TCGA Dataset Annotation context.

    Example:
        dataset=TCGA_BRCA_mRNA

    Resolves:
        TCGA_DATASET_ANNOTATIONS_DIR/
            TCGA_BRCA/
                TCGA_BRCA_CMdrug_result/
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


class TCGADatasetAnnotationCMScoreOptionsView(
    TCGADatasetAnnotationCMScoreMixin,
    BaseDatasetAnnotationCMScoreOptionsView,
):
    pass


class TCGADatasetAnnotationCMScoreResultView(
    TCGADatasetAnnotationCMScoreMixin,
    BaseDatasetAnnotationCMScoreResultView,
):
    pass


class TIMEDBDatasetAnnotationCMScoreMixin:
    """
    TIMEDB Dataset Annotation context.

    Current project group_type values:
        other | grade | stage

    The group_by/group_type pair is validated against the same
    candidate model used by Dataset Annotation availability.
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

    def get_group_by(
        self,
        request,
    ) -> str:
        group_by = str(
            request.query_params.get(
                "group_by",
                "",
            )
            or ""
        ).strip()

        if not group_by:
            raise DatasetAnnotationInputError(
                "Missing required parameter: group_by."
            )

        return group_by

    def get_group_type(
        self,
        request,
    ) -> str:
        return get_timedb_group_type_query(
            request
        )

    def validate_group_context(
        self,
        *,
        dataset_name: str,
        group_by: str | None,
        group_type: str | None,
    ) -> tuple[str, str]:
        matched = validate_timedb_group_selection(
            dataset_name=dataset_name,
            group_by=group_by,
            group_type=group_type,
        )

        return (
            matched["value"],
            matched["group_type"],
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


class TIMEDBDatasetAnnotationCMScoreOptionsView(
    TIMEDBDatasetAnnotationCMScoreMixin,
    BaseDatasetAnnotationCMScoreOptionsView,
):
    pass


class TIMEDBDatasetAnnotationCMScoreResultView(
    TIMEDBDatasetAnnotationCMScoreMixin,
    BaseDatasetAnnotationCMScoreResultView,
):
    pass


class SCSTDatasetAnnotationCMScoreMixin:
    """
    SC/ST Dataset Annotation context.

    Query params:
        dataset
        data_type
        group_by
        group_value

    Directory:
        {dataset}_{normalized_group_by}/
            {dataset}_CMdrug_result_{group_value}/
    """

    source = "SCST"

    network_source_task_type = (
        "SCSTHybridReferenceTask"
    )

    uses_generic_annotation_context = False

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

        cm_results_dir = (
            get_scst_dataset_cmdrug_result_dir_path(
                annotation_dir=annotation_dir,
                dataset_name=dataset_name,
                group_value=group_value,
            )
        )

        return {
            "dataset_name": (
                dataset_name
            ),

            "data_type": (
                data_type
            ),

            "group_by": (
                group_by
            ),

            "group_type": None,

            "group_value": (
                group_value
            ),

            "annotation_dir_name": (
                annotation_dir.name
            ),

            "annotation_file_prefix": (
                dataset_name
            ),

            "annotation_dir": (
                annotation_dir
            ),

            "cm_results_dir": (
                cm_results_dir
            ),
        }


class SCSTDatasetAnnotationCMScoreOptionsView(
    SCSTDatasetAnnotationCMScoreMixin,
    BaseDatasetAnnotationCMScoreOptionsView,
):
    pass


class SCSTDatasetAnnotationCMScoreResultView(
    SCSTDatasetAnnotationCMScoreMixin,
    BaseDatasetAnnotationCMScoreResultView,
):
    pass
