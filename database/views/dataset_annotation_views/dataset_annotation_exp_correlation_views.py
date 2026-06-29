import traceback
from pathlib import Path

import pandas as pd

from django.conf import settings
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from database.models import DatasetMetadata

from analysis.utils.paired_cohort_task_utils import (
    PAIRED_COHORT_EXPR_SAMPLE_COL,
    PAIRED_COHORT_META_SAMPLE_COL,
    PAIRED_COHORT_GROUP_COL,
    PAIRED_COHORT_CASE_LABEL,
)
from analysis.utils.hybrid_reference_task_utils import (
    HYBRID_REFERENCE_VALID_TCGA_TYPES,
)
from analysis.utils.workflow_detail_utils.workflow_exp_correlation_utils import (
    PAIRED_COHORT_VALID_EXP_CORRELATION_TYPES,
    HYBRID_REFERENCE_VALID_EXP_CORRELATION_TYPES,
    WORKFLOW_EXPR_SAMPLE_COL,
    WorkflowExpCorrelationPathError,
    WorkflowExpCorrelationInputError,
    read_exp_correlation_file_by_path,
    get_workflow_available_exp_correlation_pairs,
    get_workflow_available_exp_correlation_types,
    get_selected_workflow_exp_correlation_pair_df,
    extract_workflow_exp_correlation_stats,
    get_workflow_exp_rna_types_by_type,
    build_expression_pair_points_from_files, filter_workflow_exp_correlation_pairs_by_expression_genes,
    read_expression_gene_set_from_file, get_workflow_exp_rna_types_for_types,
    get_workflow_available_exp_correlation_types_from_pairs,
)
from database.utils.dataset_annotation_utils.path_utils import (
    DatasetAnnotationInputError,
    DatasetAnnotationPathError,
    get_dataset_query_name,
    resolve_dataset_annotation_dir,
    resolve_tcga_annotation_dir_name,
    resolve_timedb_annotation_dir_name,
    get_dataset_annotation_exp_correlation_file_path,
    validate_annotation_dataset_name,
)


RNA_TYPE_EXPRESSION_FILENAME_MAP = {
    "mRNA": "{tcga_type}_mRNA_log2tpm_exp.csv",
    "miRNA": "{tcga_type}_miRNA_log2rpm_exp.csv",
    "lncRNA": "{tcga_type}_lncRNA_log2tpm_exp.csv",
    "circRNA": "{tcga_type}_circRNA_count_exp.csv",
}

TCGA_META_FILENAME_TEMPLATE = "{tcga_type}_meta.csv"


def normalize_tcga_type(value: str) -> str:
    value = str(value or "").strip().upper()

    if not value:
        raise DatasetAnnotationInputError(
            "Missing cancer_type for TCGA expression lookup."
        )

    if value.startswith("TCGA_"):
        tcga_type = value
    else:
        tcga_type = f"TCGA_{value}"

    validate_annotation_dataset_name(tcga_type)

    if tcga_type not in HYBRID_REFERENCE_VALID_TCGA_TYPES:
        raise DatasetAnnotationInputError(
            f"Unsupported TCGA cancer type: {tcga_type}."
        )

    return tcga_type


def infer_tcga_type_from_annotation_dir_name(annotation_dir_name: str) -> str:
    """
    TCGA_ACC -> TCGA_ACC
    TCGA_ACC_mRNA -> TCGA_ACC, if it appears accidentally
    """

    annotation_dir_name = str(annotation_dir_name or "").strip()

    if annotation_dir_name.startswith("TCGA_"):
        parts = annotation_dir_name.split("_")

        if len(parts) >= 2:
            return normalize_tcga_type("_".join(parts[:2]))

    return normalize_tcga_type(annotation_dir_name)


def get_dataset_metadata_for_annotation(dataset_name: str):
    """
    TIMEDB dataset usually uses dataset name like GSE19750.
    TCGA may use TCGA_ACC_mRNA; fallback handles both.
    """

    candidates = [
        dataset_name,
        f"{dataset_name}_mRNA",
    ]

    metadata = (
        DatasetMetadata.objects
        .filter(dataset__in=candidates)
        .first()
    )

    if metadata is not None:
        return metadata

    return (
        DatasetMetadata.objects
        .filter(dataset__startswith=f"{dataset_name}_")
        .first()
    )


def get_tcga_reference_expression_file_path(
    *,
    tcga_type: str,
    rna_type: str,
) -> Path:
    rna_type = str(rna_type or "").strip()

    if rna_type not in RNA_TYPE_EXPRESSION_FILENAME_MAP:
        raise DatasetAnnotationInputError(
            "Invalid RNA type for TCGA expression lookup. "
            "Allowed values are: "
            f"{', '.join(RNA_TYPE_EXPRESSION_FILENAME_MAP.keys())}."
        )

    base_dir = Path(settings.TCGA_DATASET_BASE_DIR).resolve()

    if not base_dir.exists() or not base_dir.is_dir():
        raise DatasetAnnotationPathError(
            "TCGA dataset base directory is not available."
        )

    filename = RNA_TYPE_EXPRESSION_FILENAME_MAP[rna_type].format(
        tcga_type=tcga_type,
    )

    file_path = (base_dir / filename).resolve()

    try:
        file_path.relative_to(base_dir)
    except ValueError:
        raise DatasetAnnotationPathError(
            "Invalid TCGA reference expression file path."
        )

    return file_path


def get_tcga_reference_meta_file_path(
    *,
    tcga_type: str,
) -> Path:
    base_dir = Path(settings.TCGA_DATASET_BASE_DIR).resolve()

    if not base_dir.exists() or not base_dir.is_dir():
        raise DatasetAnnotationPathError(
            "TCGA dataset base directory is not available."
        )

    filename = TCGA_META_FILENAME_TEMPLATE.format(
        tcga_type=tcga_type,
    )

    file_path = (base_dir / filename).resolve()

    try:
        file_path.relative_to(base_dir)
    except ValueError:
        raise DatasetAnnotationPathError(
            "Invalid TCGA reference meta file path."
        )

    return file_path


class BaseDatasetAnnotationExpCorrelationMixin:
    source = None
    network_source_task_type = None

    annotation_root_setting_name = None
    annotation_dir_name_resolver = None

    valid_types = []

    def get_expression_file_path(
            self,
            *,
            context: dict,
            rna_type: str,
    ) -> Path:
        return get_tcga_reference_expression_file_path(
            tcga_type=context["tcga_type"],
            rna_type=rna_type,
        )

    def get_expression_gene_sets(
            self,
            context: dict,
    ) -> dict[str, set[str]]:
        expression_gene_sets = {}

        for rna_type in get_workflow_exp_rna_types_for_types(self.valid_types):
            try:
                expression_file = self.get_expression_file_path(
                    context=context,
                    rna_type=rna_type,
                )

                expression_gene_sets[rna_type] = (
                    read_expression_gene_set_from_file(
                        expression_file,
                        sample_col=WORKFLOW_EXPR_SAMPLE_COL,
                    )
                )

            except FileNotFoundError:
                expression_gene_sets[rna_type] = set()

        return expression_gene_sets

    def filter_pairs_by_expression_genes(
            self,
            *,
            context: dict,
            pairs: list[dict],
    ) -> list[dict]:
        expression_gene_sets = self.get_expression_gene_sets(context)

        return filter_workflow_exp_correlation_pairs_by_expression_genes(
            pairs=pairs,
            expression_gene_sets=expression_gene_sets,
        )

    def get_annotation_root_dir(self):
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

        correlation_file = get_dataset_annotation_exp_correlation_file_path(
            annotation_dir=annotation_dir,
            file_prefix=file_prefix,
        )

        metadata = get_dataset_metadata_for_annotation(dataset_name)

        context = {
            "dataset_name": dataset_name,
            "annotation_dir_name": annotation_dir_name,
            "annotation_file_prefix": file_prefix,
            "annotation_dir": annotation_dir,
            "correlation_file": correlation_file,
            "metadata": metadata,
        }

        context["tcga_type"] = self.get_tcga_type(context)

        return context

    def get_tcga_type(self, context: dict) -> str:
        raise NotImplementedError

    def get_base_response(self, context: dict) -> dict:
        return {
            "success": True,
            "source": self.source,
            "dataset_name": context["dataset_name"],
            "annotation_dir_name": context["annotation_dir_name"],
            "annotation_file_prefix": context["annotation_file_prefix"],
            "network_source_task_type": self.network_source_task_type,
            "tcga_type": context["tcga_type"],
        }


class DatasetAnnotationExpCorrelationOptionsBaseView(
    BaseDatasetAnnotationExpCorrelationMixin,
    APIView,
):
    """
    Return selectable gene1/gene2/type pairs for dataset annotation
    expression correlation plot.

    Query params:
        dataset or dataset_name

    Input file:
        {annotation_file_prefix}_ceRNA_corr.csv
    """

    def get(self, request):
        try:
            if not self.source:
                raise RuntimeError("Missing source.")

            if not self.network_source_task_type:
                raise RuntimeError("Missing network_source_task_type.")

            if not self.valid_types:
                raise RuntimeError("Missing valid_types.")

            context = self.resolve_annotation_context(request)

            try:
                correlation_file, cor_df = read_exp_correlation_file_by_path(
                    context["correlation_file"]
                )

                raw_results = get_workflow_available_exp_correlation_pairs(
                    df=cor_df,
                    valid_types=self.valid_types,
                )

                results = self.filter_pairs_by_expression_genes(
                    context=context,
                    pairs=raw_results,
                )

                available_types = get_workflow_available_exp_correlation_types_from_pairs(
                    pairs=results,
                    valid_types=self.valid_types,
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
                WorkflowExpCorrelationInputError,
                WorkflowExpCorrelationPathError,
            ) as e:
                return Response(
                    {
                        "success": False,
                        "detail": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            response_data = {
                **self.get_base_response(context),
                "correlation_file": correlation_file.name,
                "valid_types": self.valid_types,
                "available_types": available_types,
                "raw_count": len(raw_results),
                "count": len(results),
                "dropped_count": len(raw_results) - len(results),
                "results": results,
            }

            return Response(
                response_data,
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


class DatasetAnnotationExpCorrelationPlotDataBaseView(
    BaseDatasetAnnotationExpCorrelationMixin,
    APIView,
):
    """
    Return scatter plot data for dataset annotation expression correlation.

    Query params:
        dataset or dataset_name
        gene1
        gene2
        type
    """

    def get_sample_ids(
        self,
        *,
        context: dict,
    ) -> list[str] | None:
        return None

    def get(self, request):
        try:
            if not self.source:
                raise RuntimeError("Missing source.")

            if not self.network_source_task_type:
                raise RuntimeError("Missing network_source_task_type.")

            if not self.valid_types:
                raise RuntimeError("Missing valid_types.")

            context = self.resolve_annotation_context(request)

            gene1 = str(
                request.query_params.get("gene1", "")
            ).strip()

            gene2 = str(
                request.query_params.get("gene2", "")
            ).strip()

            type_value = str(
                request.query_params.get("type", "")
            ).strip()

            if not gene1:
                return Response(
                    {
                        "success": False,
                        "detail": "Missing query parameter: gene1.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not gene2:
                return Response(
                    {
                        "success": False,
                        "detail": "Missing query parameter: gene2.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not type_value:
                return Response(
                    {
                        "success": False,
                        "detail": "Missing query parameter: type.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if type_value not in self.valid_types:
                return Response(
                    {
                        "success": False,
                        "detail": (
                            "Invalid type. Allowed values are: "
                            f"{', '.join(self.valid_types)}."
                        ),
                        "valid_types": self.valid_types,
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                correlation_file, cor_df = read_exp_correlation_file_by_path(
                    context["correlation_file"]
                )

                pair_df = get_selected_workflow_exp_correlation_pair_df(
                    cor_df=cor_df,
                    gene1=gene1,
                    gene2=gene2,
                    type_value=type_value,
                )

                correlation_stats = extract_workflow_exp_correlation_stats(
                    pair_df
                )

                rna_type_map = get_workflow_exp_rna_types_by_type(type_value)

                gene1_expr_file = self.get_expression_file_path(
                    context=context,
                    rna_type=rna_type_map["gene1_rna_type"],
                )

                gene2_expr_file = self.get_expression_file_path(
                    context=context,
                    rna_type=rna_type_map["gene2_rna_type"],
                )

                sample_ids = self.get_sample_ids(
                    context=context,
                )

                point_result = build_expression_pair_points_from_files(
                    gene1_expr_file=gene1_expr_file,
                    gene2_expr_file=gene2_expr_file,
                    gene1=gene1,
                    gene2=gene2,
                    sample_col=WORKFLOW_EXPR_SAMPLE_COL,
                    sample_ids=sample_ids,
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
                WorkflowExpCorrelationInputError,
                WorkflowExpCorrelationPathError,
                DatasetAnnotationInputError,
                DatasetAnnotationPathError,
                ValueError,
            ) as e:
                return Response(
                    {
                        "success": False,
                        "detail": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            response_data = {
                **self.get_base_response(context),
                "type": type_value,
                "gene1": gene1,
                "gene2": gene2,
                "correlation_file": correlation_file.name,
                "gene1_expression_file": gene1_expr_file.name,
                "gene2_expression_file": gene2_expr_file.name,
                "correlation": correlation_stats,
            }

            response_data.update(point_result)

            return Response(
                response_data,
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


class TCGADatasetAnnotationExpCorrelationMixin:
    source = "TCGA"
    network_source_task_type = "PairedCohortTask"

    annotation_root_setting_name = "TCGA_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_tcga_annotation_dir_name
    )

    valid_types = PAIRED_COHORT_VALID_EXP_CORRELATION_TYPES

    def get_tcga_type(self, context: dict) -> str:
        return infer_tcga_type_from_annotation_dir_name(
            context["annotation_dir_name"]
        )


class TCGADatasetAnnotationExpCorrelationOptionsView(
    TCGADatasetAnnotationExpCorrelationMixin,
    DatasetAnnotationExpCorrelationOptionsBaseView,
):
    pass


class TCGADatasetAnnotationExpCorrelationPlotDataView(
    TCGADatasetAnnotationExpCorrelationMixin,
    DatasetAnnotationExpCorrelationPlotDataBaseView,
):
    def get_sample_ids(
        self,
        *,
        context: dict,
    ) -> list[str] | None:
        meta_file = get_tcga_reference_meta_file_path(
            tcga_type=context["tcga_type"]
        )

        if not meta_file.exists() or not meta_file.is_file():
            raise FileNotFoundError(
                f"TCGA meta file not found: {meta_file.name}."
            )

        try:
            meta_df = pd.read_csv(
                meta_file,
                index_col=PAIRED_COHORT_META_SAMPLE_COL,
            )
        except ValueError as e:
            raise WorkflowExpCorrelationInputError(
                f"Failed to read TCGA meta file by sample column "
                f"{PAIRED_COHORT_META_SAMPLE_COL}: {str(e)}"
            )
        except Exception as e:
            raise WorkflowExpCorrelationInputError(
                f"Failed to read TCGA meta file: {str(e)}"
            )

        if PAIRED_COHORT_GROUP_COL not in meta_df.columns:
            raise WorkflowExpCorrelationInputError(
                f"TCGA meta file is missing required column: "
                f"{PAIRED_COHORT_GROUP_COL}."
            )

        case_samples = (
            meta_df[
                meta_df[PAIRED_COHORT_GROUP_COL]
                .astype(str)
                .str.strip()
                == PAIRED_COHORT_CASE_LABEL
            ]
            .index.astype(str)
            .tolist()
        )

        if not case_samples:
            raise WorkflowExpCorrelationInputError(
                f"No case samples found in TCGA meta file by "
                f"{PAIRED_COHORT_GROUP_COL}={PAIRED_COHORT_CASE_LABEL}."
            )

        return case_samples


class TIMEDBDatasetAnnotationExpCorrelationMixin:
    source = "TIMEDB"
    network_source_task_type = "HybridReferenceTask"

    annotation_root_setting_name = "TIMEDB_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_timedb_annotation_dir_name
    )

    valid_types = HYBRID_REFERENCE_VALID_EXP_CORRELATION_TYPES

    def get_tcga_type(self, context: dict) -> str:
        metadata = context.get("metadata")

        if metadata is None:
            raise DatasetAnnotationInputError(
                "Dataset metadata is required to resolve TIMEDB cancer_type."
            )

        return normalize_tcga_type(
            getattr(metadata, "cancer_type", "")
        )


class TIMEDBDatasetAnnotationExpCorrelationOptionsView(
    TIMEDBDatasetAnnotationExpCorrelationMixin,
    DatasetAnnotationExpCorrelationOptionsBaseView,
):
    pass


class TIMEDBDatasetAnnotationExpCorrelationPlotDataView(
    TIMEDBDatasetAnnotationExpCorrelationMixin,
    DatasetAnnotationExpCorrelationPlotDataBaseView,
):
    pass
