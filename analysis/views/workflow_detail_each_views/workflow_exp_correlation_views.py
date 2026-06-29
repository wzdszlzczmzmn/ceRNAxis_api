import traceback

import pandas as pd

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.models import (
    PairedCohortTask,
    HybridReferenceTask,
)
from analysis.utils.hybrid_reference_task_utils import get_hybrid_reference_tcga_expression_file_path

from analysis.utils.workflow_detail_utils.workflow_network_view_utils import (
    WorkflowNetworkViewError,
    get_required_task_uuid,
    get_task_or_raise,
    require_success_task,
)

from analysis.utils.workflow_detail_utils.workflow_exp_correlation_utils import (
    PAIRED_COHORT_VALID_EXP_CORRELATION_TYPES,
    HYBRID_REFERENCE_VALID_EXP_CORRELATION_TYPES,
    WorkflowExpCorrelationPathError,
    WorkflowExpCorrelationInputError,
    read_workflow_exp_correlation_file,
    get_workflow_available_exp_correlation_pairs,
    get_workflow_available_exp_correlation_types,
    get_workflow_exp_file_fields_by_type,
    get_selected_workflow_exp_correlation_pair_df,
    extract_workflow_exp_correlation_stats,
    get_workflow_exp_rna_types_by_type,
    build_expression_pair_points_from_files, get_workflow_exp_rna_types_for_types, read_expression_gene_set_from_file,
    filter_workflow_exp_correlation_pairs_by_expression_genes, get_workflow_available_exp_correlation_types_from_pairs,
)

from analysis.utils.paired_cohort_task_utils import (
    PAIRED_COHORT_EXPR_SAMPLE_COL,
    PAIRED_COHORT_META_SAMPLE_COL,
    PAIRED_COHORT_GROUP_COL,
    PAIRED_COHORT_CASE_LABEL,
    PairedCohortTaskPathError, get_paired_cohort_input_file_path,
)


class WorkflowExpCorrelationOptionsBaseView(APIView):
    """
    Return all selectable gene1/gene2/type pairs for workflow expression
    correlation plot.

    Query params:
        taskUUID: Workflow task UUID

    Input file:
        {task_name}_ceRNA_corr.csv
    """

    task_model = None
    task_type = None
    task_label = "Workflow task"
    valid_types = []

    should_filter_pairs_by_expression_genes = True

    def get_extra_response_data(self, task) -> dict:
        return {}

    def get_expression_file_path_by_rna_type(
        self,
        *,
        task,
        rna_type: str,
    ):
        raise NotImplementedError

    def get_expression_gene_sets(self, task) -> dict[str, set[str]]:
        expression_gene_sets = {}

        for rna_type in get_workflow_exp_rna_types_for_types(self.valid_types):
            try:
                expression_file = self.get_expression_file_path_by_rna_type(
                    task=task,
                    rna_type=rna_type,
                )

                expression_gene_sets[rna_type] = (
                    read_expression_gene_set_from_file(
                        expression_file,
                    )
                )

            except FileNotFoundError:
                expression_gene_sets[rna_type] = set()

        return expression_gene_sets

    def filter_pairs(
        self,
        *,
        task,
        pairs: list[dict],
    ) -> list[dict]:
        if not self.should_filter_pairs_by_expression_genes:
            return pairs

        expression_gene_sets = self.get_expression_gene_sets(task)

        return filter_workflow_exp_correlation_pairs_by_expression_genes(
            pairs=pairs,
            expression_gene_sets=expression_gene_sets,
        )

    def get(self, request):
        try:
            if self.task_model is None:
                raise RuntimeError("Missing task_model.")

            if not self.task_type:
                raise RuntimeError("Missing task_type.")

            if not self.valid_types:
                raise RuntimeError("Missing valid_types.")

            task_uuid = get_required_task_uuid(request)

            task = get_task_or_raise(
                model_class=self.task_model,
                task_uuid=task_uuid,
                task_label=self.task_type,
            )

            require_success_task(
                task=task,
                task_label=self.task_label,
            )

            try:
                correlation_file, cor_df = read_workflow_exp_correlation_file(
                    task
                )

                raw_results = get_workflow_available_exp_correlation_pairs(
                    df=cor_df,
                    valid_types=self.valid_types,
                )

                results = self.filter_pairs(
                    task=task,
                    pairs=raw_results,
                )

                available_types = get_workflow_available_exp_correlation_types_from_pairs(
                    pairs=results,
                    valid_types=self.valid_types,
                )

            except FileNotFoundError as e:
                return Response(
                    {
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
                        "detail": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            response_data = {
                "uuid": str(task.uuid),
                "task_type": self.task_type,
                "task_name": task.task_name,
                "correlation_file": correlation_file.name,
                "valid_types": self.valid_types,
                "available_types": available_types,
                "raw_count": len(raw_results),
                "count": len(results),
                "dropped_count": len(raw_results) - len(results),
                "results": results,
            }

            response_data.update(
                self.get_extra_response_data(task)
            )

            return Response(
                response_data,
                status=status.HTTP_200_OK,
            )

        except WorkflowNetworkViewError as e:
            return Response(
                {
                    "detail": e.msg,
                },
                status=e.status_code,
            )

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {
                    "detail": f"Server error: {str(e)}",
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )


class WorkflowExpCorrelationPlotDataBaseView(APIView):
    """
    Return scatter plot data for workflow expression correlation plot.

    Query params:
        taskUUID: Workflow task UUID
        gene1: selected gene1
        gene2: selected gene2
        type: selected interaction type
    """

    task_model = None
    task_type = None
    task_label = "Workflow task"
    valid_types = []

    def get_extra_response_data(self, task) -> dict:
        return {}

    def build_plot_points(
        self,
        *,
        task,
        cor_df: pd.DataFrame,
        pair_df: pd.DataFrame,
        gene1: str,
        gene2: str,
        type_value: str,
        correlation_file_name: str,
    ) -> dict:
        """
        Subclass hook.

        Must return:
            summary, regression, points
        """
        raise NotImplementedError

    def get(self, request):
        try:
            if self.task_model is None:
                raise RuntimeError("Missing task_model.")

            if not self.task_type:
                raise RuntimeError("Missing task_type.")

            if not self.valid_types:
                raise RuntimeError("Missing valid_types.")

            task_uuid = get_required_task_uuid(request)

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
                        "detail": "Missing query parameter: gene1."
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not gene2:
                return Response(
                    {
                        "detail": "Missing query parameter: gene2."
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not type_value:
                return Response(
                    {
                        "detail": "Missing query parameter: type."
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if type_value not in self.valid_types:
                return Response(
                    {
                        "detail": (
                            "Invalid type. Allowed values are: "
                            f"{', '.join(self.valid_types)}."
                        ),
                        "valid_types": self.valid_types,
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            task = get_task_or_raise(
                model_class=self.task_model,
                task_uuid=task_uuid,
                task_label=self.task_type,
            )

            require_success_task(
                task=task,
                task_label=self.task_label,
            )

            try:
                correlation_file, cor_df = read_workflow_exp_correlation_file(
                    task
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

                point_result = self.build_plot_points(
                    task=task,
                    cor_df=cor_df,
                    pair_df=pair_df,
                    gene1=gene1,
                    gene2=gene2,
                    type_value=type_value,
                    correlation_file_name=correlation_file.name,
                )

            except FileNotFoundError as e:
                return Response(
                    {
                        "detail": str(e),
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )
            except (
                WorkflowExpCorrelationInputError,
                WorkflowExpCorrelationPathError,
                PairedCohortTaskPathError,
                ValueError,
            ) as e:
                return Response(
                    {
                        "detail": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )
            except Exception as e:
                return Response(
                    {
                        "detail": f"Failed to build plot data: {str(e)}"
                    },
                    status=status.HTTP_500_INTERNAL_SERVER_ERROR,
                )

            response_data = {
                "uuid": str(task.uuid),
                "task_type": self.task_type,
                "task_name": task.task_name,
                "type": type_value,
                "gene1": gene1,
                "gene2": gene2,
                "correlation_file": correlation_file.name,
                "correlation": correlation_stats,
            }

            response_data.update(point_result)
            response_data.update(self.get_extra_response_data(task))

            return Response(
                response_data,
                status=status.HTTP_200_OK,
            )

        except WorkflowNetworkViewError as e:
            return Response(
                {
                    "detail": e.msg,
                },
                status=e.status_code,
            )

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {
                    "detail": f"Server error: {str(e)}",
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )


class PairedCohortExpCorrelationOptionsView(
    WorkflowExpCorrelationOptionsBaseView
):
    task_model = PairedCohortTask
    task_type = "PairedCohortTask"
    task_label = "Paired cohort task"
    valid_types = PAIRED_COHORT_VALID_EXP_CORRELATION_TYPES

    RNA_TYPE_FILE_FIELD_MAP = {
        "mRNA": "mrna_file",
        "miRNA": "mirna_file",
        "lncRNA": "lncrna_file",
        "circRNA": "circrna_file",
    }

    def get_expression_file_path_by_rna_type(
        self,
        *,
        task,
        rna_type: str,
    ):
        field_name = self.RNA_TYPE_FILE_FIELD_MAP.get(rna_type)

        if not field_name:
            raise WorkflowExpCorrelationInputError(
                f"Unsupported RNA type: {rna_type}."
            )

        return get_paired_cohort_input_file_path(
            task,
            field_name,
        )

    def get_extra_response_data(self, task) -> dict:
        return {
            "map_info": task.map_info,
            "deg_method": task.deg_method,
            "use_padj": getattr(task, "use_padj", True),
        }


class PairedCohortExpCorrelationPlotDataView(
    WorkflowExpCorrelationPlotDataBaseView
):
    task_model = PairedCohortTask
    task_type = "PairedCohortTask"
    task_label = "Paired cohort task"
    valid_types = PAIRED_COHORT_VALID_EXP_CORRELATION_TYPES

    def build_plot_points(
        self,
        *,
        task,
        cor_df: pd.DataFrame,
        pair_df: pd.DataFrame,
        gene1: str,
        gene2: str,
        type_value: str,
        correlation_file_name: str,
    ) -> dict:
        file_map = get_workflow_exp_file_fields_by_type(type_value)

        gene1_expr_file = get_paired_cohort_input_file_path(
            task,
            file_map["gene1_file"],
        )

        gene2_expr_file = get_paired_cohort_input_file_path(
            task,
            file_map["gene2_file"],
        )

        meta_file = get_paired_cohort_input_file_path(
            task,
            "meta_file",
        )

        if not meta_file.exists() or not meta_file.is_file():
            raise FileNotFoundError(
                f"Input file not found: {meta_file.name}."
            )

        try:
            meta_df = pd.read_csv(
                meta_file,
                index_col=PAIRED_COHORT_META_SAMPLE_COL,
            )
        except ValueError as e:
            raise WorkflowExpCorrelationInputError(
                f"Failed to read meta file by sample column "
                f"{PAIRED_COHORT_META_SAMPLE_COL}: {str(e)}"
            )
        except Exception as e:
            raise WorkflowExpCorrelationInputError(
                f"Failed to read meta file: {str(e)}"
            )

        if PAIRED_COHORT_GROUP_COL not in meta_df.columns:
            raise WorkflowExpCorrelationInputError(
                f"Meta file is missing required column: "
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
                f"No case samples found in meta file by "
                f"{PAIRED_COHORT_GROUP_COL}={PAIRED_COHORT_CASE_LABEL}."
            )

        return build_expression_pair_points_from_files(
            gene1_expr_file=gene1_expr_file,
            gene2_expr_file=gene2_expr_file,
            gene1=gene1,
            gene2=gene2,
            sample_col=PAIRED_COHORT_EXPR_SAMPLE_COL,
            sample_ids=case_samples,
        )

    def get_extra_response_data(self, task) -> dict:
        return {
            "map_info": task.map_info,
            "deg_method": task.deg_method,
            "use_padj": getattr(task, "use_padj", True),
        }

class HybridReferenceExpCorrelationOptionsView(
    WorkflowExpCorrelationOptionsBaseView
):
    task_model = HybridReferenceTask
    task_type = "HybridReferenceTask"
    task_label = "Hybrid reference task"
    valid_types = HYBRID_REFERENCE_VALID_EXP_CORRELATION_TYPES

    def get_expression_file_path_by_rna_type(
        self,
        *,
        task,
        rna_type: str,
    ):
        return get_hybrid_reference_tcga_expression_file_path(
            task=task,
            rna_type=rna_type,
        )

    def get_extra_response_data(self, task) -> dict:
        return {
            "map_info": task.map_info,
            "tcga_type": task.tcga_type,
            "lncrna_type": task.lncrna_type,
            "deg_method": task.deg_method,
            "use_padj": getattr(task, "use_padj", True),
        }


class HybridReferenceExpCorrelationPlotDataView(
    WorkflowExpCorrelationPlotDataBaseView
):
    task_model = HybridReferenceTask
    task_type = "HybridReferenceTask"
    task_label = "Hybrid reference task"
    valid_types = HYBRID_REFERENCE_VALID_EXP_CORRELATION_TYPES

    def build_plot_points(
        self,
        *,
        task,
        cor_df: pd.DataFrame,
        pair_df: pd.DataFrame,
        gene1: str,
        gene2: str,
        type_value: str,
        correlation_file_name: str,
    ) -> dict:
        rna_type_map = get_workflow_exp_rna_types_by_type(type_value)

        gene1_rna_type = rna_type_map["gene1_rna_type"]
        gene2_rna_type = rna_type_map["gene2_rna_type"]

        gene1_expr_file = get_hybrid_reference_tcga_expression_file_path(
            task=task,
            rna_type=gene1_rna_type,
        )

        gene2_expr_file = get_hybrid_reference_tcga_expression_file_path(
            task=task,
            rna_type=gene2_rna_type,
        )

        return build_expression_pair_points_from_files(
            gene1_expr_file=gene1_expr_file,
            gene2_expr_file=gene2_expr_file,
            gene1=gene1,
            gene2=gene2,
            sample_col="sample_id",
            sample_ids=None,
        )

    def get_extra_response_data(self, task) -> dict:
        return {
            "map_info": task.map_info,
            "tcga_type": task.tcga_type,
            "lncrna_type": task.lncrna_type,
            "deg_method": task.deg_method,
            "use_padj": getattr(task, "use_padj", True),
        }
