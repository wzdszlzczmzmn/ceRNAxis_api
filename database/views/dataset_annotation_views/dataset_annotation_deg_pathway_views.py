import traceback

import pandas as pd

from django.conf import settings

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.utils.workflow_detail_utils.deg_pathway_utils import (
    DEGPathwayInputError,
    validate_deg_pathway_dataframe_columns,
    build_deg_pathway_data_from_dataframe_common,
)
from database.utils.dataset_annotation_utils.path_utils import (
    DatasetAnnotationInputError,
    DatasetAnnotationPathError,
    get_dataset_query_name,
    resolve_dataset_annotation_dir,
    resolve_tcga_annotation_dir_name,
    resolve_timedb_annotation_dir_name,
    get_dataset_annotation_mrna_gsea_file_path, get_timedb_group_type_query, resolve_timedb_group_annotation_dir_name,
    resolve_timedb_group_annotation_file_prefix,
)


class BaseDatasetAnnotationDEGPathwayView(APIView):
    """
    Return DEG pathway bubble plot data for dataset annotation.

    Query params:
        dataset or dataset_name

    Input file:
        {annotation_file_prefix}_mRNA_gsea.csv
    """

    source = None
    network_source_task_type = None

    annotation_root_setting_name = None
    annotation_dir_name_resolver = None

    title = "DEG Pathway Enrichment"

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

        gsea_file = get_dataset_annotation_mrna_gsea_file_path(
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
            "gsea_file": gsea_file,
        }

    def read_gsea_file(self, gsea_file):
        if not gsea_file.exists() or not gsea_file.is_file():
            raise FileNotFoundError(
                f"Dataset annotation mRNA GSEA file not found: "
                f"{gsea_file.name}"
            )

        try:
            df = pd.read_csv(gsea_file)
        except Exception as e:
            raise DEGPathwayInputError(
                f"Failed to read dataset annotation mRNA GSEA file: "
                f"{str(e)}"
            )

        validate_deg_pathway_dataframe_columns(df=df)

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

            if not self.annotation_dir_name_resolver:
                raise RuntimeError("Missing annotation_dir_name_resolver.")

            context = self.resolve_annotation_context(request)

            try:
                df = self.read_gsea_file(
                    context["gsea_file"]
                )

                result = build_deg_pathway_data_from_dataframe_common(
                    gsea_file_name=context["gsea_file"].name,
                    df=df,
                    title=self.title,
                    base_response=self.get_base_response(context),
                )

            except FileNotFoundError as e:
                return Response(
                    {
                        "success": False,
                        "detail": str(e),
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            except DEGPathwayInputError as e:
                return Response(
                    {
                        "success": False,
                        "detail": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            return Response(
                result,
                status=status.HTTP_200_OK,
            )

        except DatasetAnnotationInputError as e:
            return Response(
                {
                    "success": False,
                    "detail": str(e),
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        except DatasetAnnotationPathError as e:
            return Response(
                {
                    "success": False,
                    "detail": str(e),
                },
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


class TCGADatasetAnnotationDEGPathwayView(
    BaseDatasetAnnotationDEGPathwayView
):
    """
    TCGA dataset annotation DEG pathway data.

    Source semantics:
        Paired Cohort annotation output.

    Input file:
        {annotation_file_prefix}_mRNA_gsea.csv
    """

    source = "TCGA"
    network_source_task_type = "PairedCohortTask"

    annotation_root_setting_name = "TCGA_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_tcga_annotation_dir_name
    )

    title = "DEG Pathway Enrichment"


class TIMEDBDatasetAnnotationDEGPathwayView(
    BaseDatasetAnnotationDEGPathwayView
):
    """
    TIMEDB dataset annotation DEG pathway data.

    Source semantics:
        Hybrid Reference annotation output.

    Input file:
        {annotation_file_prefix}_mRNA_gsea.csv
    """

    source = "TIMEDB"
    network_source_task_type = "HybridReferenceTask"

    annotation_root_setting_name = "TIMEDB_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_timedb_annotation_dir_name
    )

    title = "TCGA-based DEG Pathway Enrichment"

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
