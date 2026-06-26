import csv
import traceback
import uuid as uuid_lib
from pathlib import Path
import pandas as pd
import numpy as np

from django.db.models import Q
from django.http import FileResponse

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.services.task_common.registry import TaskNotFoundError, MultipleTaskMatchedError
from analysis.services.task_download.registry import TaskDownloadConfigNotFoundError
from analysis.services.task_download.service import prepare_task_result_download, InvalidTaskUUIDError, \
    TaskNotReadyForDownloadError, TaskResultFileNotFoundError, TaskResultArchiveError, TaskDownloadError
from analysis.utils.paired_cohort_network_utils import paired_cohort_rna_file_node_key, \
    PAIRED_COHORT_IMMUNE_AXIS_REQUIRED_COLUMNS, parse_database_list, paired_cohort_cerna_edge_id, \
    append_database_to_edge, PAIRED_COHORT_CERNA_TYPE_TO_TARGET_RNA_TYPE, PAIRED_COHORT_VALID_CERNA_EDGE_TYPES, \
    PAIRED_COHORT_CERNA_AXIS_REQUIRED_COLUMNS, get_paired_cohort_immune_axis_file_path, \
    get_paired_cohort_cerna_axis_file_path
from analysis.utils.paired_cohort_task_utils import PAIRED_COHORT_INPUT_FILENAME_MAP, \
    get_paired_cohort_input_file_path, PairedCohortTaskPathError, get_paired_cohort_task_output_dir, validate_safe_name, \
    read_paired_cohort_correlation_file, get_paired_cohort_available_correlation_pairs, PairedCohortTaskInputError, \
    PAIRED_COHORT_VALID_CORRELATION_TYPES, normalize_paired_cohort_correlation_dataframe, \
    get_paired_cohort_exp_file_fields_by_type, PAIRED_COHORT_EXPR_SAMPLE_COL, PAIRED_COHORT_META_SAMPLE_COL, \
    PAIRED_COHORT_GROUP_COL, PAIRED_COHORT_CASE_LABEL, extract_paired_cohort_correlation_stats, \
    read_paired_cohort_background_file, get_paired_cohort_available_background_types, \
    read_paired_cohort_axis_final_file, normalize_paired_cohort_axis_final_dataframe, PAIRED_COHORT_AXIS_FINAL_COLUMNS, \
    build_paired_cohort_survival_km_data
from analysis.utils.workflow_network_utils import immune_gene_node_key, normalize_rna_name, \
    append_immune_annotation_to_edge, immune_edge_id, edge_pair_key, format_immune_evidence_item, \
    mark_node_as_immune_related, to_float_or_none, is_valid_uuid
from database.models import RNANode, RNAInteraction
from analysis.models import CustomListQueryTask, PairedCohortTask

from analysis.utils.custom_list_query_task_utils import (
    get_task_output_dir,
    get_task_out_prefix,
)


class CustomListQueryTaskNetworkView(APIView):
    """
    根据 CustomListQueryTask 中保存的 RNA 列表查询 workflow network。

    GET 参数：
        taskUUID=<uuid>

    网络组成：
    1. 用户提交 RNA 列表中已存在于 RNANode 的节点
    2. 这些节点之间的 RNAInteraction 边
    3. output/{task_name}_map_immune_gene_axis.csv 中的 immune annotation 信息

    immune 合并规则：
    1. immune source 已由脚本保证一定来自 task.rnas["miRNA"]
    2. immune target 如果已经是当前网络中的节点，则复用并标记 is_immune_related
    3. immune target 如果不在当前网络中，则新增 immune_checkpoint_gene 节点
    4. immune edge 如果与已有 RNAInteraction edge 重叠，则给已有边追加 immune_annotation
    5. immune edge 如果不与已有边重叠，则新增 immune_annotation edge
    """

    ALLOWED_RNA_TYPES = [
        "miRNA",
        "mRNA",
        "lncRNA",
        "circRNA",
    ]

    IMMUNE_RESULT_REQUIRED_COLUMNS = {
        "miRNA",
        "Immune checkpointGene",
    }

    def get(self, request, *args, **kwargs):
        try:
            task_uuid = request.query_params.get("taskUUID", "").strip()

            if not task_uuid:
                return Response(
                    {
                        "success": False,
                        "msg": "Missing required parameter: taskUUID.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not self.is_valid_uuid(task_uuid):
                return Response(
                    {
                        "success": False,
                        "msg": "Invalid taskUUID format.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                task = CustomListQueryTask.objects.get(uuid=task_uuid)
            except CustomListQueryTask.DoesNotExist:
                return Response(
                    {
                        "success": False,
                        "msg": f"CustomListQueryTask with UUID {task_uuid} not found.",
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            input_items = self.parse_task_rnas(task.rnas)

            if not input_items:
                return Response(
                    self.empty_network_response(task),
                    status=status.HTTP_200_OK,
                )

            result = self.build_network(
                task=task,
                input_items=input_items,
            )

            result["task_type"] = "CustomListQueryTask"
            result["task_uuid"] = str(task.uuid)
            result["task_name"] = task.task_name
            result["map_info"] = task.map_info

            return Response(result, status=status.HTTP_200_OK)

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {
                    "success": False,
                    "msg": f"Server error: {str(e)}",
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

    @staticmethod
    def is_valid_uuid(value: str) -> bool:
        try:
            uuid_lib.UUID(str(value))
            return True
        except ValueError:
            return False

    def empty_network_response(self, task):
        return {
            "task_type": "CustomListQueryTask",
            "task_uuid": str(task.uuid),
            "task_name": task.task_name,
            "map_info": task.map_info,
            "meta": {
                "input_node_count": 0,
                "matched_node_count": 0,
                "returned_node_count": 0,
                "ignored_node_count": 0,
                "edge_count": 0,
                "rna_interaction_edge_count": 0,
                "immune_annotation_edge_count": 0,
                "immune_merged_edge_count": 0,
                "immune_new_edge_count": 0,
                "immune_annotation_loaded": False,
                "edge_scope": "matched_input_nodes_only_with_immune_annotation",
            },
            "nodes": [],
            "edges": [],
            "ignored_nodes": [],
        }

    def parse_task_rnas(self, rnas):
        """
        将 task.rnas 转成：
        [
            {"name": "MALAT1", "type": "lncRNA"},
            ...
        ]
        """
        if not isinstance(rnas, dict):
            return []

        items = []
        seen = set()

        for rna_type in self.ALLOWED_RNA_TYPES:
            names = rnas.get(rna_type, [])

            if not isinstance(names, list):
                continue

            for name in names:
                if not isinstance(name, str):
                    continue

                name = name.strip()

                if not name:
                    continue

                key = (name, rna_type)

                if key in seen:
                    continue

                seen.add(key)

                items.append(
                    {
                        "name": name,
                        "type": rna_type,
                    }
                )

        return items

    def build_network(self, task, input_items):
        input_key_set = {
            (item["name"], item["type"])
            for item in input_items
        }

        input_mirna_name_map = {
            self.normalize_rna_name(item["name"]): item["name"]
            for item in input_items
            if item["type"] == "miRNA"
        }

        matched_nodes = self.query_matched_nodes(input_key_set)

        matched_node_map = {
            (node.name, node.rna_type): node
            for node in matched_nodes
        }

        matched_node_name_map = {}

        for node in matched_nodes:
            matched_node_name_map.setdefault(node.name, []).append(node)

        matched_node_ids = [
            node.id
            for node in matched_nodes
        ]

        ignored_nodes = self.build_ignored_nodes_for_unmatched_inputs(
            input_items=input_items,
            matched_node_map=matched_node_map,
        )

        nodes = {}
        edges = {}
        edge_pair_index = {}
        used_node_keys = set()

        rna_interaction_edge_count = self.add_rna_interaction_edges(
            matched_node_ids=matched_node_ids,
            nodes=nodes,
            edges=edges,
            edge_pair_index=edge_pair_index,
            used_node_keys=used_node_keys,
        )

        immune_result = self.add_immune_annotation_edges(
            task=task,
            input_mirna_name_map=input_mirna_name_map,
            nodes=nodes,
            edges=edges,
            edge_pair_index=edge_pair_index,
            used_node_keys=used_node_keys,
            matched_node_map=matched_node_map,
            matched_node_name_map=matched_node_name_map,
        )

        for node in matched_nodes:
            node_key = self.rna_node_key(node)

            if node_key not in used_node_keys:
                ignored_nodes.append(
                    {
                        "name": node.name,
                        "type": node.rna_type,
                        "reason": "no_interaction_in_query_scope",
                    }
                )

        return {
            "meta": {
                "input_node_count": len(input_items),
                "matched_node_count": len(matched_nodes),
                "returned_node_count": len(nodes),
                "ignored_node_count": len(ignored_nodes),
                "edge_count": len(edges),
                "rna_interaction_edge_count": rna_interaction_edge_count,
                "immune_annotation_edge_count": immune_result["edge_count"],
                "immune_merged_edge_count": immune_result["merged_edge_count"],
                "immune_new_edge_count": immune_result["new_edge_count"],
                "immune_annotation_loaded": immune_result["loaded"],
                "immune_annotation_file": immune_result["file"],
                "edge_scope": "matched_input_nodes_only_with_immune_annotation",
            },
            "nodes": list(nodes.values()),
            "edges": list(edges.values()),
            "ignored_nodes": ignored_nodes,
        }

    def query_matched_nodes(self, input_key_set):
        if not input_key_set:
            return []

        node_query = Q()

        for name, rna_type in input_key_set:
            node_query |= Q(
                name=name,
                rna_type=rna_type,
            )

        return list(
            RNANode.objects.filter(node_query)
        )

    def build_ignored_nodes_for_unmatched_inputs(
            self,
            input_items,
            matched_node_map,
    ):
        ignored_nodes = []

        for item in input_items:
            key = (item["name"], item["type"])

            if key not in matched_node_map:
                ignored_nodes.append(
                    {
                        "name": item["name"],
                        "type": item["type"],
                        "reason": "not_found_in_database",
                    }
                )

        return ignored_nodes

    def add_rna_interaction_edges(
            self,
            matched_node_ids,
            nodes,
            edges,
            edge_pair_index,
            used_node_keys,
    ) -> int:
        if matched_node_ids:
            interactions = (
                RNAInteraction.objects
                .filter(
                    source_id__in=matched_node_ids,
                    target_id__in=matched_node_ids,
                )
                .select_related("source", "target")
                .prefetch_related("databases")
                .order_by("id")
            )
        else:
            interactions = RNAInteraction.objects.none()

        edge_count = 0

        for interaction in interactions:
            source = interaction.source
            target = interaction.target

            source_key = self.rna_node_key(source)
            target_key = self.rna_node_key(target)

            used_node_keys.add(source_key)
            used_node_keys.add(target_key)

            if source_key not in nodes:
                nodes[source_key] = self.format_rna_node(source)

            if target_key not in nodes:
                nodes[target_key] = self.format_rna_node(target)

            edge_id = f"rna_interaction:{interaction.id}"

            if edge_id in edges:
                continue

            edges[edge_id] = {
                "id": edge_id,
                "source": source_key,
                "target": target_key,
                "source_db_id": source.id,
                "target_db_id": target.id,
                "source_name": source.name,
                "target_name": target.name,
                "source_type": source.rna_type,
                "target_type": target.rna_type,
                "species": interaction.species,
                "type": interaction.interaction_type,
                "edge_type": "rna_interaction",
                "databases": [
                    db.name
                    for db in interaction.databases.all()
                ],
                "is_immune_related": False,
                "immune_annotation": None,
            }

            edge_pair_key = self.edge_pair_key(source_key, target_key)

            if edge_pair_key not in edge_pair_index:
                edge_pair_index[edge_pair_key] = edge_id

            edge_count += 1

        return edge_count

    def add_immune_annotation_edges(
            self,
            task,
            input_mirna_name_map,
            nodes,
            edges,
            edge_pair_index,
            used_node_keys,
            matched_node_map,
            matched_node_name_map,
    ) -> dict:
        immune_file_path = self.get_immune_result_file_path(task)

        if not immune_file_path.exists() or not immune_file_path.is_file():
            return {
                "loaded": False,
                "edge_count": 0,
                "merged_edge_count": 0,
                "new_edge_count": 0,
                "file": str(immune_file_path),
            }

        total_immune_records = 0
        merged_edge_count = 0
        new_edge_count = 0

        with immune_file_path.open("r", encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f)

            if not reader.fieldnames:
                return {
                    "loaded": False,
                    "edge_count": 0,
                    "merged_edge_count": 0,
                    "new_edge_count": 0,
                    "file": str(immune_file_path),
                }

            missing_columns = (
                    self.IMMUNE_RESULT_REQUIRED_COLUMNS
                    - set(reader.fieldnames)
            )

            if missing_columns:
                raise ValueError(
                    f"Missing required immune result columns: {sorted(missing_columns)}"
                )

            for row in reader:
                mirna_name = str(row.get("miRNA", "")).strip()
                immune_gene = str(row.get("Immune checkpointGene", "")).strip()

                if not mirna_name or not immune_gene:
                    continue

                normalized_mirna_name = self.normalize_rna_name(mirna_name)

                if normalized_mirna_name not in input_mirna_name_map:
                    continue

                canonical_mirna_name = input_mirna_name_map[normalized_mirna_name]

                source_key = self.get_or_create_mirna_node(
                    mirna_name=canonical_mirna_name,
                    nodes=nodes,
                    matched_node_map=matched_node_map,
                )

                target_key = self.get_or_create_immune_target_node(
                    immune_gene=immune_gene,
                    nodes=nodes,
                    matched_node_map=matched_node_map,
                    matched_node_name_map=matched_node_name_map,
                )

                self.mark_node_as_immune_related(
                    nodes=nodes,
                    node_key=source_key,
                    immune_node_type="immune_source_miRNA",
                    pathway=row.get("Immune checkpointPathway"),
                )

                self.mark_node_as_immune_related(
                    nodes=nodes,
                    node_key=target_key,
                    immune_node_type="immune_checkpoint_gene",
                    pathway=row.get("Immune checkpointPathway"),
                )

                used_node_keys.add(source_key)
                used_node_keys.add(target_key)

                evidence_item = self.format_immune_evidence_item(row)

                edge_pair_key = self.edge_pair_key(source_key, target_key)

                if edge_pair_key in edge_pair_index:
                    existing_edge_id = edge_pair_index[edge_pair_key]
                    self.append_immune_annotation_to_edge(
                        edge=edges[existing_edge_id],
                        evidence_item=evidence_item,
                    )
                    merged_edge_count += 1
                else:
                    edge_id = self.immune_edge_id(
                        mirna_name=canonical_mirna_name,
                        immune_gene=immune_gene,
                    )

                    if edge_id not in edges:
                        edges[edge_id] = {
                            "id": edge_id,
                            "source": source_key,
                            "target": target_key,
                            "source_name": canonical_mirna_name,
                            "target_name": immune_gene,
                            "source_type": "miRNA",
                            "target_type": "immune_checkpoint_gene",
                            "species": "Homo sapiens",
                            "type": "immune_annotation",
                            "edge_type": "immune_annotation",
                            "databases": ["immune_annotation"],
                            "is_immune_related": True,
                            "immune_annotation": {
                                "evidence_count": 0,
                                "evidence_items": [],
                            },
                        }

                        edge_pair_index[edge_pair_key] = edge_id
                        new_edge_count += 1

                    self.append_immune_annotation_to_edge(
                        edge=edges[edge_id],
                        evidence_item=evidence_item,
                    )

                total_immune_records += 1

        return {
            "loaded": True,
            "edge_count": total_immune_records,
            "merged_edge_count": merged_edge_count,
            "new_edge_count": new_edge_count,
            "file": str(immune_file_path),
        }

    def get_or_create_mirna_node(
            self,
            mirna_name,
            nodes,
            matched_node_map,
    ) -> str:
        matched_node = matched_node_map.get(
            (mirna_name, "miRNA")
        )

        if matched_node is not None:
            node_key = self.rna_node_key(matched_node)

            if node_key not in nodes:
                nodes[node_key] = self.format_rna_node(matched_node)

            return node_key

        node_key = self.task_mirna_node_key(mirna_name)

        if node_key not in nodes:
            nodes[node_key] = {
                "id": node_key,
                "db_id": None,
                "name": mirna_name,
                "type": "miRNA",
                "species": None,
                "source": "task_input",
                "matched_in_database": False,
                "is_immune_related": True,
                "immune_node_type": "immune_source_miRNA",
                "immune_pathways": [],
            }

        return node_key

    def get_or_create_immune_target_node(
            self,
            immune_gene,
            nodes,
            matched_node_map,
            matched_node_name_map,
    ) -> str:
        """
        优先复用用户提交 RNA 中已经匹配到数据库的同名节点。

        例如：
        - 用户提交 mRNA 里有 CD274，并且 RNANode 中存在 CD274/mRNA
        - immune target 也是 CD274
        则复用 node:<db_id>，不新建 immune_checkpoint_gene:CD274。
        """

        preferred_types = [
            "mRNA",
            "lncRNA",
            "circRNA",
            "miRNA",
        ]

        for rna_type in preferred_types:
            matched_node = matched_node_map.get(
                (immune_gene, rna_type)
            )

            if matched_node is not None:
                node_key = self.rna_node_key(matched_node)

                if node_key not in nodes:
                    nodes[node_key] = self.format_rna_node(matched_node)

                return node_key

        same_name_nodes = matched_node_name_map.get(immune_gene, [])

        if same_name_nodes:
            matched_node = same_name_nodes[0]
            node_key = self.rna_node_key(matched_node)

            if node_key not in nodes:
                nodes[node_key] = self.format_rna_node(matched_node)

            return node_key

        node_key = self.immune_gene_node_key(immune_gene)

        if node_key not in nodes:
            nodes[node_key] = {
                "id": node_key,
                "db_id": None,
                "name": immune_gene,
                "type": "immune_checkpoint_gene",
                "species": "Homo sapiens",
                "source": "immune_annotation",
                "matched_in_database": False,
                "is_immune_related": True,
                "immune_node_type": "immune_checkpoint_gene",
                "immune_pathways": [],
            }

        return node_key

    @staticmethod
    def mark_node_as_immune_related(
            nodes,
            node_key,
            immune_node_type,
            pathway=None,
    ):
        if node_key not in nodes:
            return

        node = nodes[node_key]

        node["is_immune_related"] = True
        node["immune_node_type"] = immune_node_type

        if "immune_pathways" not in node or not isinstance(
                node["immune_pathways"],
                list,
        ):
            node["immune_pathways"] = []

        pathway = str(pathway or "").strip()

        if pathway and pathway not in node["immune_pathways"]:
            node["immune_pathways"].append(pathway)

    def append_immune_annotation_to_edge(
            self,
            edge,
            evidence_item,
    ):
        edge["is_immune_related"] = True

        if edge.get("immune_annotation") is None:
            edge["immune_annotation"] = {
                "evidence_count": 0,
                "evidence_items": [],
            }

        evidence_items = edge["immune_annotation"]["evidence_items"]

        evidence_key = self.immune_evidence_key(evidence_item)

        existing_keys = {
            self.immune_evidence_key(item)
            for item in evidence_items
        }

        if evidence_key not in existing_keys:
            evidence_items.append(evidence_item)

        edge["immune_annotation"]["evidence_count"] = len(evidence_items)

    def format_immune_evidence_item(self, row):
        return {
            "correlation_coefficient": self.to_float_or_none(
                row.get("CorrelationCoefficient")
            ),
            "p_value": self.to_float_or_none(
                row.get("P-value")
            ),
            "q_value": self.to_float_or_none(
                row.get("Q-value")
            ),
            "cancer": str(row.get("Cancer", "")).strip(),
            "pathway": str(row.get("Immune checkpointPathway", "")).strip(),
            "evidence": str(row.get("Evidence", "")).strip(),
        }

    @staticmethod
    def immune_evidence_key(evidence_item):
        return (
            evidence_item.get("correlation_coefficient"),
            evidence_item.get("p_value"),
            evidence_item.get("q_value"),
            evidence_item.get("cancer"),
            evidence_item.get("pathway"),
            evidence_item.get("evidence"),
        )

    def get_immune_result_file_path(self, task) -> Path:
        output_dir = get_task_output_dir(task)
        out_prefix = get_task_out_prefix(task)

        return output_dir / f"{out_prefix}.csv"

    @staticmethod
    def to_float_or_none(value):
        if value is None:
            return None

        value = str(value).strip()

        if not value:
            return None

        try:
            return float(value)
        except ValueError:
            return None

    @staticmethod
    def format_rna_node(node):
        node_key = f"node:{node.id}"

        return {
            "id": node_key,
            "db_id": node.id,
            "name": node.name,
            "type": node.rna_type,
            "species": node.species,
            "source": "rna_interaction_database",
            "matched_in_database": True,
            "is_immune_related": False,
            "immune_node_type": None,
            "immune_pathways": [],
        }

    @staticmethod
    def rna_node_key(node):
        return f"node:{node.id}"

    @staticmethod
    def task_mirna_node_key(mirna_name: str):
        return f"task_miRNA:{mirna_name}"

    @staticmethod
    def immune_gene_node_key(gene_name: str):
        return f"immune_checkpoint_gene:{gene_name}"

    @staticmethod
    def edge_pair_key(source_key: str, target_key: str):
        return f"{source_key}=>{target_key}"

    @staticmethod
    def immune_edge_id(mirna_name: str, immune_gene: str):
        return f"immune_annotation:{mirna_name}:{immune_gene}"

    @staticmethod
    def normalize_rna_name(name: str) -> str:
        return str(name or "").strip().lower()


class PairedCohortUploadedFileDownloadView(APIView):
    """
    Download uploaded input file for PairedCohortTask.

    Query params:
        taskUUID: PairedCohortTask UUID
        file_type: one of mrna_file, mirna_file, lncrna_file, meta_file
    """

    def get(self, request):
        try:
            task_uuid = str(request.query_params.get("taskUUID", "")).strip()
            file_type = str(request.query_params.get("file_type", "")).strip()

            if not task_uuid:
                return Response(
                    {
                        "success": False,
                        "msg": "Missing query parameter: taskUUID.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                uuid_lib.UUID(task_uuid)
            except ValueError:
                return Response(
                    {
                        "success": False,
                        "msg": "Invalid Task UUID format.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not file_type:
                return Response(
                    {
                        "success": False,
                        "msg": "Missing query parameter: file_type.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if file_type not in PAIRED_COHORT_INPUT_FILENAME_MAP:
                return Response(
                    {
                        "success": False,
                        "msg": (
                            "Invalid file_type. Allowed values are: "
                            f"{', '.join(PAIRED_COHORT_INPUT_FILENAME_MAP.keys())}."
                        ),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                task = PairedCohortTask.objects.get(uuid=task_uuid)
            except PairedCohortTask.DoesNotExist:
                return Response(
                    {
                        "success": False,
                        "msg": f"PairedCohortTask not found: {task_uuid}.",
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            try:
                file_path = get_paired_cohort_input_file_path(
                    task=task,
                    field_name=file_type,
                )
            except PairedCohortTaskPathError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not file_path.exists() or not file_path.is_file():
                return Response(
                    {
                        "success": False,
                        "msg": (
                            "Uploaded file not found: "
                            f"{PAIRED_COHORT_INPUT_FILENAME_MAP[file_type]}."
                        ),
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            return FileResponse(
                open(file_path, "rb"),
                as_attachment=True,
                filename=file_path.name,
                content_type="text/csv",
            )

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {
                    "success": False,
                    "msg": f"Server error: {str(e)}",
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )


class PairedCohortDEGVolcanoView(APIView):
    """
    Return DEG volcano plot data for PairedCohortTask.

    Query params:
        taskUUID: PairedCohortTask UUID
        rna_type: one of mRNA, miRNA, lncRNA, circRNA

    DEG filename rule:
        {task_name}_{deg_method}_{rna_type}.csv

    Required DEG columns:
        gene_name, log2FC, pvalue, regulation
    """

    REQUIRED_COLUMNS = {
        "gene_name",
        "log2FC",
        "pvalue",
        "regulation",
    }

    VALID_RNA_TYPES = ["mRNA", "miRNA", "lncRNA", "circRNA"]

    VALID_REGULATION_GROUPS = ["NotSig", "Down", "Up"]

    def get(self, request):
        try:
            task_uuid = str(request.query_params.get("taskUUID", "")).strip()
            rna_type = str(request.query_params.get("rna_type", "")).strip()

            if not task_uuid:
                return Response(
                    {"detail": "Missing query parameter: taskUUID."},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                uuid_lib.UUID(task_uuid)
            except ValueError:
                return Response(
                    {"detail": "Invalid Task UUID format."},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not rna_type:
                return Response(
                    {"detail": "Missing query parameter: rna_type."},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if rna_type not in self.VALID_RNA_TYPES:
                return Response(
                    {
                        "detail": (
                            "Invalid rna_type. Allowed values are: "
                            f"{', '.join(self.VALID_RNA_TYPES)}."
                        )
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                task = PairedCohortTask.objects.get(uuid=task_uuid)
            except PairedCohortTask.DoesNotExist:
                return Response(
                    {
                        "detail": f"PairedCohortTask not found: {task_uuid}."
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            if task.status != PairedCohortTask.Status.Success:
                return Response(
                    {
                        "detail": (
                            "Paired cohort task is not completed successfully. "
                            f"Current status: {task.get_status_display()}."
                        )
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                deg_file = self.get_deg_file_path(
                    task=task,
                    rna_type=rna_type,
                )
            except PairedCohortTaskPathError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not deg_file.exists() or not deg_file.is_file():
                return Response(
                    {
                        "detail": f"DEG file not found: {deg_file.name}."
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            try:
                df = pd.read_csv(deg_file)
            except Exception as e:
                return Response(
                    {"detail": f"Failed to read DEG file: {e}"},
                    status=status.HTTP_500_INTERNAL_SERVER_ERROR,
                )

            missing_columns = self.REQUIRED_COLUMNS - set(df.columns)

            if missing_columns:
                return Response(
                    {
                        "detail": (
                            "Missing required columns: "
                            f"{sorted(missing_columns)}"
                        )
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            df = df[
                [
                    "gene_name",
                    "log2FC",
                    "pvalue",
                    "regulation",
                ]
            ].copy()

            df = df.replace([np.inf, -np.inf], np.nan)

            raw_count = int(df.shape[0])

            df = df.dropna(
                subset=[
                    "gene_name",
                    "log2FC",
                    "pvalue",
                    "regulation",
                ]
            )

            df["log2FC"] = pd.to_numeric(
                df["log2FC"],
                errors="coerce",
            )

            df["pvalue"] = pd.to_numeric(
                df["pvalue"],
                errors="coerce",
            )

            df = df.dropna(
                subset=[
                    "log2FC",
                    "pvalue",
                ]
            )

            df = df[df["pvalue"] > 0]

            df = df[
                df["regulation"].isin(self.VALID_REGULATION_GROUPS)
            ]

            cleaned_count = int(df.shape[0])
            dropped_count = raw_count - cleaned_count

            df["neg_log10_pvalue"] = -np.log10(df["pvalue"])

            groups = {}

            for group in self.VALID_REGULATION_GROUPS:
                sub_df = df[df["regulation"] == group]

                groups[group] = [
                    {
                        "gene_name": row["gene_name"],
                        "log2FC": float(row["log2FC"]),
                        "pvalue": float(row["pvalue"]),
                        "neg_log10_pvalue": float(row["neg_log10_pvalue"]),
                    }
                    for _, row in sub_df.iterrows()
                ]

            return Response(
                {
                    "uuid": str(task.uuid),
                    "task_name": task.task_name,
                    "deg_method": task.deg_method,
                    "rna_type": rna_type,
                    "deg_file": deg_file.name,
                    "summary": {
                        "raw_count": raw_count,
                        "cleaned_count": cleaned_count,
                        "dropped_count": dropped_count,
                        "not_sig": len(groups["NotSig"]),
                        "down": len(groups["Down"]),
                        "up": len(groups["Up"]),
                    },
                    "groups": groups,
                },
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

    @staticmethod
    def get_deg_file_path(task, rna_type: str):
        task_name = str(task.task_name).strip()
        deg_method = str(task.deg_method).strip()

        validate_safe_name(task_name, "task_name")
        validate_safe_name(deg_method, "deg_method")

        output_dir = get_paired_cohort_task_output_dir(task)

        filename = f"{task_name}_{deg_method}_{rna_type}.csv"
        deg_file = (output_dir / filename).resolve()

        if not str(deg_file).startswith(str(output_dir)):
            raise PairedCohortTaskPathError(
                "Invalid paired cohort DEG file path."
            )

        return deg_file


class PairedCohortLog2FCCorrelationView(APIView):
    """
        Return core data for log2FC correlation plot.

        Query params:
            taskUUID: PairedCohortTask UUID
            type: one of miRNA-mRNA, miRNA-lncRNA

        Input filename:
            {task_name}_ceRNA_background.csv
        """

    REQUIRED_COLUMNS = {
        "miRNA",
        "ceRNA",
        "species",
        "database",
        "type",
        "miRNA_log2FC",
        "ceRNA_log2FC",
    }

    X_COL = "ceRNA_log2FC"
    Y_COL = "miRNA_log2FC"

    def get(self, request):
        try:
            task_uuid = str(request.query_params.get("taskUUID", "")).strip()
            requested_type = request.query_params.get("type", None)

            type_value = (
                str(requested_type).strip()
                if requested_type is not None
                else ""
            )

            if not task_uuid:
                return Response(
                    {"detail": "Missing query parameter: taskUUID."},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                uuid_lib.UUID(task_uuid)
            except ValueError:
                return Response(
                    {"detail": "Invalid Task UUID format."},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                task = PairedCohortTask.objects.get(uuid=task_uuid)
            except PairedCohortTask.DoesNotExist:
                return Response(
                    {
                        "detail": f"PairedCohortTask not found: {task_uuid}."
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            if task.status != PairedCohortTask.Status.Success:
                return Response(
                    {
                        "detail": (
                            "Paired cohort task is not completed successfully. "
                            f"Current status: {task.get_status_display()}."
                        )
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                background_file, df = read_paired_cohort_background_file(task)
            except PairedCohortTaskPathError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_400_BAD_REQUEST,
                )
            except FileNotFoundError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_404_NOT_FOUND,
                )
            except PairedCohortTaskInputError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            available_types = get_paired_cohort_available_background_types(df)

            if not available_types:
                return Response(
                    {
                        "detail": (
                            "No supported background interaction type found "
                            "in ceRNA background file."
                        )
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not type_value:
                type_value = available_types[0]

            if type_value not in available_types:
                return Response(
                    {
                        "detail": (
                            "Invalid type for this task. Allowed values are: "
                            f"{', '.join(available_types)}."
                        ),
                        "available_types": available_types,
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            result = self.build_response_data(
                task=task,
                df=df,
                type_value=type_value,
                background_file_name=background_file.name,
                available_types=available_types,
            )

            return Response(result, status=status.HTTP_200_OK)

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {"detail": f"Server error: {str(e)}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )


    def build_response_data(
            self,
            task,
            df: pd.DataFrame,
            type_value: str,
            background_file_name: str,
            available_types: list[str],
    ) -> dict:
        x_col = self.X_COL
        y_col = self.Y_COL

        df = df[df["type"] == type_value].copy()

        raw_count = int(df.shape[0])

        df[x_col] = pd.to_numeric(df[x_col], errors="coerce")
        df[y_col] = pd.to_numeric(df[y_col], errors="coerce")

        required_subset = [
            "miRNA",
            "ceRNA",
            "type",
            x_col,
            y_col,
        ]

        df = df.replace([np.inf, -np.inf], np.nan)
        df = df.dropna(subset=required_subset).copy()

        cleaned_count = int(df.shape[0])
        dropped_count = raw_count - cleaned_count

        if cleaned_count == 0:
            return {
                "uuid": str(task.uuid),
                "task_name": task.task_name,
                "type": type_value,
                "available_types": available_types,
                "background_file": background_file_name,
                "summary": {
                    "raw_count": raw_count,
                    "cleaned_count": 0,
                    "dropped_count": dropped_count,
                    "anti_count": 0,
                    "same_count": 0,
                },
                "points": [],
            }

        df["anti_correlation"] = df[x_col] * df[y_col] < 0

        anti_count = int(df["anti_correlation"].sum())
        same_count = cleaned_count - anti_count

        return {
            "uuid": str(task.uuid),
            "task_name": task.task_name,
            "type": type_value,
            "available_types": available_types,
            "background_file": background_file_name,
            "summary": {
                "raw_count": raw_count,
                "cleaned_count": cleaned_count,
                "dropped_count": dropped_count,
                "anti_count": anti_count,
                "same_count": same_count,
            },
            "points": self.serialize_points(df),
        }

    def serialize_points(self, df: pd.DataFrame) -> list[dict]:
        points = []

        for _, row in df.iterrows():
            species = row.get("species", "")
            database = row.get("database", "")

            if pd.isna(species):
                species = ""

            if pd.isna(database):
                database = ""

            points.append(
                {
                    "miRNA": row["miRNA"],
                    "ceRNA": row["ceRNA"],
                    "species": species,
                    "database": database,
                    "type": row["type"],
                    "ceRNA_log2FC": float(row[self.X_COL]),
                    "miRNA_log2FC": float(row[self.Y_COL]),
                    "anti_correlation": bool(row["anti_correlation"]),
                }
            )

        return points


class PairedCohortExpCorrelationOptionsView(APIView):
    """
    Return all selectable gene1/gene2/type pairs for expression correlation plot.

    Query params:
        taskUUID: PairedCohortTask UUID
    """

    def get(self, request):
        try:
            task_uuid = str(request.query_params.get("taskUUID", "")).strip()

            if not task_uuid:
                return Response(
                    {"detail": "Missing query parameter: taskUUID."},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                uuid_lib.UUID(task_uuid)
            except ValueError:
                return Response(
                    {"detail": "Invalid Task UUID format."},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                task = PairedCohortTask.objects.get(uuid=task_uuid)
            except PairedCohortTask.DoesNotExist:
                return Response(
                    {"detail": f"PairedCohortTask not found: {task_uuid}."},
                    status=status.HTTP_404_NOT_FOUND,
                )

            if task.status != PairedCohortTask.Status.Success:
                return Response(
                    {
                        "detail": (
                            "Paired cohort task is not completed successfully. "
                            f"Current status: {task.get_status_display()}."
                        )
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                correlation_file, cor_df = read_paired_cohort_correlation_file(task)
                results = get_paired_cohort_available_correlation_pairs(cor_df)
            except FileNotFoundError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_404_NOT_FOUND,
                )
            except (PairedCohortTaskInputError, PairedCohortTaskPathError) as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            return Response(
                {
                    "uuid": str(task.uuid),
                    "task_name": task.task_name,
                    "correlation_file": correlation_file.name,
                    "valid_types": PAIRED_COHORT_VALID_CORRELATION_TYPES,
                    "count": len(results),
                    "results": results,
                },
                status=status.HTTP_200_OK,
            )

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {"detail": f"Server error: {str(e)}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )


class PairedCohortExpCorrelationPlotDataView(APIView):
    """
    Return core scatter plot data for expression correlation plot.

    Query params:
        taskUUID: PairedCohortTask UUID
        gene1: selected gene1
        gene2: selected gene2
        type: one of miRNA-mRNA, miRNA-lncRNA, lncRNA-mRNA
    """

    def get(self, request):
        try:
            task_uuid = str(request.query_params.get("taskUUID", "")).strip()
            gene1 = str(request.query_params.get("gene1", "")).strip()
            gene2 = str(request.query_params.get("gene2", "")).strip()
            type_value = str(request.query_params.get("type", "")).strip()

            if not task_uuid:
                return Response(
                    {"detail": "Missing query parameter: taskUUID."},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not gene1:
                return Response(
                    {"detail": "Missing query parameter: gene1."},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not gene2:
                return Response(
                    {"detail": "Missing query parameter: gene2."},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not type_value:
                return Response(
                    {"detail": "Missing query parameter: type."},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                uuid_lib.UUID(task_uuid)
            except ValueError:
                return Response(
                    {"detail": "Invalid Task UUID format."},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if type_value not in PAIRED_COHORT_VALID_CORRELATION_TYPES:
                return Response(
                    {
                        "detail": (
                            "Invalid type. Allowed values are: "
                            f"{', '.join(PAIRED_COHORT_VALID_CORRELATION_TYPES)}."
                        )
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                task = PairedCohortTask.objects.get(uuid=task_uuid)
            except PairedCohortTask.DoesNotExist:
                return Response(
                    {"detail": f"PairedCohortTask not found: {task_uuid}."},
                    status=status.HTTP_404_NOT_FOUND,
                )

            if task.status != PairedCohortTask.Status.Success:
                return Response(
                    {
                        "detail": (
                            "Paired cohort task is not completed successfully. "
                            f"Current status: {task.get_status_display()}."
                        )
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                correlation_file, cor_df = read_paired_cohort_correlation_file(task)

                result = self.build_plot_data(
                    task=task,
                    cor_df=cor_df,
                    gene1=gene1,
                    gene2=gene2,
                    type_value=type_value,
                    correlation_file_name=correlation_file.name,
                )

            except FileNotFoundError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_404_NOT_FOUND,
                )
            except (PairedCohortTaskInputError, PairedCohortTaskPathError, ValueError) as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_400_BAD_REQUEST,
                )
            except Exception as e:
                return Response(
                    {"detail": f"Failed to build plot data: {str(e)}"},
                    status=status.HTTP_500_INTERNAL_SERVER_ERROR,
                )

            return Response(result, status=status.HTTP_200_OK)

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {"detail": f"Server error: {str(e)}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

    @staticmethod
    def build_plot_data(
            task,
            cor_df: pd.DataFrame,
            gene1: str,
            gene2: str,
            type_value: str,
            correlation_file_name: str,
    ) -> dict:
        normalized_cor_df = normalize_paired_cohort_correlation_dataframe(cor_df)

        pair_df = normalized_cor_df[
            (normalized_cor_df["gene1"] == gene1)
            & (normalized_cor_df["gene2"] == gene2)
            & (normalized_cor_df["type"] == type_value)
            ].copy()

        if pair_df.empty:
            raise ValueError(
                "Selected gene pair was not found in paired cohort correlation file."
            )

        file_map = get_paired_cohort_exp_file_fields_by_type(type_value)

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

        for file_path in [gene1_expr_file, gene2_expr_file, meta_file]:
            if not file_path.exists() or not file_path.is_file():
                raise FileNotFoundError(
                    f"Input file not found: {file_path.name}."
                )

        gene1_df = pd.read_csv(
            gene1_expr_file,
            index_col=PAIRED_COHORT_EXPR_SAMPLE_COL,
        )
        gene2_df = pd.read_csv(
            gene2_expr_file,
            index_col=PAIRED_COHORT_EXPR_SAMPLE_COL,
        )
        meta_df = pd.read_csv(
            meta_file,
            index_col=PAIRED_COHORT_META_SAMPLE_COL,
        )

        if gene1 not in gene1_df.columns:
            raise ValueError(
                f"gene1 not found in expression file: {gene1}."
            )

        if gene2 not in gene2_df.columns:
            raise ValueError(
                f"gene2 not found in expression file: {gene2}."
            )

        if PAIRED_COHORT_GROUP_COL not in meta_df.columns:
            raise ValueError(
                f"Meta file is missing required column: {PAIRED_COHORT_GROUP_COL}."
            )

        case_samples = (
            meta_df[
                meta_df[PAIRED_COHORT_GROUP_COL].astype(str).str.strip()
                == PAIRED_COHORT_CASE_LABEL
                ]
            .index.astype(str)
            .tolist()
        )

        if not case_samples:
            raise ValueError(
                f"No case samples found in meta file by "
                f"{PAIRED_COHORT_GROUP_COL}={PAIRED_COHORT_CASE_LABEL}."
            )

        gene1_expr = gene1_df[[gene1]].copy()
        gene2_expr = gene2_df[[gene2]].copy()

        gene1_expr.index = gene1_expr.index.astype(str)
        gene2_expr.index = gene2_expr.index.astype(str)

        merged_df = pd.merge(
            gene1_expr,
            gene2_expr,
            left_index=True,
            right_index=True,
            how="inner",
        )

        merged_df = merged_df[merged_df.index.isin(case_samples)].copy()
        merged_df.columns = ["gene1_expr", "gene2_expr"]

        merged_df["gene1_expr"] = pd.to_numeric(
            merged_df["gene1_expr"],
            errors="coerce",
        )
        merged_df["gene2_expr"] = pd.to_numeric(
            merged_df["gene2_expr"],
            errors="coerce",
        )

        raw_count = int(merged_df.shape[0])

        plot_df = (
            merged_df
            .replace([np.inf, -np.inf], np.nan)
            .dropna(subset=["gene1_expr", "gene2_expr"])
            .copy()
        )

        cleaned_count = int(plot_df.shape[0])
        dropped_count = raw_count - cleaned_count

        correlation_stats = extract_paired_cohort_correlation_stats(pair_df)
        regression = PairedCohortExpCorrelationPlotDataView.calculate_regression(
            plot_df
        )

        points = [
            {
                "sample_id": str(sample_id),
                "gene1_expr": float(row["gene1_expr"]),
                "gene2_expr": float(row["gene2_expr"]),
            }
            for sample_id, row in plot_df.iterrows()
        ]

        return {
            "uuid": str(task.uuid),
            "task_name": task.task_name,
            "type": type_value,
            "gene1": gene1,
            "gene2": gene2,
            "correlation_file": correlation_file_name,
            "summary": {
                "raw_count": raw_count,
                "cleaned_count": cleaned_count,
                "dropped_count": dropped_count,
            },
            "correlation": correlation_stats,
            "regression": regression,
            "points": points,
        }

    @staticmethod
    def calculate_regression(plot_df: pd.DataFrame) -> dict:
        if plot_df.shape[0] < 2:
            return {
                "slope": None,
                "intercept": None,
            }

        if plot_df["gene1_expr"].nunique() < 2:
            return {
                "slope": None,
                "intercept": None,
            }

        slope, intercept = np.polyfit(
            plot_df["gene1_expr"],
            plot_df["gene2_expr"],
            1,
        )

        return {
            "slope": float(slope),
            "intercept": float(intercept),
        }


class PairedCohortTaskNetworkView(APIView):
    """
    Build network for PairedCohortTask from workflow output files.

    GET params:
        taskUUID=<uuid>

    Network source files:
        output/{task_name}_ceRNA_axis.csv
        output/{task_name}_map_immune_axis.csv

    Response contract:
        The node / edge structure is aligned with CustomListQueryTaskNetworkView
        so that the same frontend network component can render both task types.
    """

    def get(self, request, *args, **kwargs):
        try:
            task_uuid = str(request.query_params.get("taskUUID", "")).strip()

            if not task_uuid:
                return Response(
                    {
                        "success": False,
                        "msg": "Missing required parameter: taskUUID.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            if not is_valid_uuid(task_uuid):
                return Response(
                    {
                        "success": False,
                        "msg": "Invalid taskUUID format.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                task = PairedCohortTask.objects.get(uuid=task_uuid)
            except PairedCohortTask.DoesNotExist:
                return Response(
                    {
                        "success": False,
                        "msg": (
                            f"PairedCohortTask with UUID {task_uuid} "
                            "not found."
                        ),
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            if task.status != PairedCohortTask.Status.Success:
                return Response(
                    {
                        "success": False,
                        "msg": (
                            "Paired cohort task is not completed successfully. "
                            f"Current status: {task.get_status_display()}."
                        ),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                result = self.build_network(task)
            except PairedCohortTaskPathError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )
            except ValueError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            result["task_type"] = "PairedCohortTask"
            result["task_uuid"] = str(task.uuid)
            result["task_name"] = task.task_name
            result["map_info"] = task.map_info
            result["deg_method"] = task.deg_method

            return Response(result, status=status.HTTP_200_OK)

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {
                    "success": False,
                    "msg": f"Server error: {str(e)}",
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

    def build_network(self, task) -> dict:
        ceRNA_axis_file = get_paired_cohort_cerna_axis_file_path(task)
        immune_axis_file = get_paired_cohort_immune_axis_file_path(task)

        if not ceRNA_axis_file.exists() or not ceRNA_axis_file.is_file():
            return self.empty_network_response(
                task=task,
                ceRNA_axis_file=ceRNA_axis_file,
                immune_axis_file=immune_axis_file,
            )

        nodes = {}
        edges = {}
        edge_pair_index = {}
        node_name_index = {}

        cerna_result = self.add_cerna_axis_edges(
            ceRNA_axis_file=ceRNA_axis_file,
            nodes=nodes,
            edges=edges,
            edge_pair_index=edge_pair_index,
            node_name_index=node_name_index,
        )

        immune_result = self.add_immune_annotation_edges(
            immune_axis_file=immune_axis_file,
            nodes=nodes,
            edges=edges,
            edge_pair_index=edge_pair_index,
            node_name_index=node_name_index,
        )

        return {
            "meta": {
                "input_node_count": None,
                "matched_node_count": None,
                "returned_node_count": len(nodes),
                "ignored_node_count": 0,
                "edge_count": len(edges),

                # 与 CustomListQueryTaskNetworkView 对齐：
                # 前端可以继续读取 rna_interaction_edge_count。
                "rna_interaction_edge_count": cerna_result["edge_count"],

                # Paired Cohort 专属统计，保留给详情页或调试使用。
                "cerna_axis_edge_count": cerna_result["edge_count"],
                "cerna_axis_loaded": cerna_result["loaded"],

                "immune_annotation_edge_count": immune_result["edge_count"],
                "immune_merged_edge_count": immune_result["merged_edge_count"],
                "immune_new_edge_count": immune_result["new_edge_count"],
                "immune_annotation_loaded": immune_result["loaded"],

                # 如前端不需要文件路径，可保留但不展示。
                "cerna_axis_file": cerna_result["file"],
                "immune_annotation_file": immune_result["file"],

                "edge_scope": "paired_cohort_output_files",
            },
            "nodes": list(nodes.values()),
            "edges": list(edges.values()),
            "ignored_nodes": [],
        }

    def empty_network_response(
            self,
            task,
            ceRNA_axis_file: Path,
            immune_axis_file: Path,
    ) -> dict:
        return {
            "task_type": "PairedCohortTask",
            "task_uuid": str(task.uuid),
            "task_name": task.task_name,
            "map_info": task.map_info,
            "deg_method": task.deg_method,
            "meta": {
                "input_node_count": None,
                "matched_node_count": None,
                "returned_node_count": 0,
                "ignored_node_count": 0,
                "edge_count": 0,
                "rna_interaction_edge_count": 0,
                "cerna_axis_edge_count": 0,
                "cerna_axis_loaded": False,
                "cerna_axis_file": str(ceRNA_axis_file),
                "immune_annotation_edge_count": 0,
                "immune_merged_edge_count": 0,
                "immune_new_edge_count": 0,
                "immune_annotation_loaded": False,
                "immune_annotation_file": str(immune_axis_file),
                "edge_scope": "paired_cohort_output_files",
            },
            "nodes": [],
            "edges": [],
            "ignored_nodes": [],
        }

    def add_cerna_axis_edges(
            self,
            ceRNA_axis_file: Path,
            nodes: dict,
            edges: dict,
            edge_pair_index: dict,
            node_name_index: dict,
    ) -> dict:
        edge_count = 0

        with ceRNA_axis_file.open("r", encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f)

            if not reader.fieldnames:
                return {
                    "loaded": False,
                    "edge_count": 0,
                    "file": str(ceRNA_axis_file),
                }

            missing_columns = (
                    PAIRED_COHORT_CERNA_AXIS_REQUIRED_COLUMNS
                    - set(reader.fieldnames)
            )

            if missing_columns:
                raise ValueError(
                    "Missing required ceRNA axis columns: "
                    f"{sorted(missing_columns)}"
                )

            for row in reader:
                mirna_name = str(row.get("miRNA", "")).strip()
                cerna_name = str(row.get("ceRNA", "")).strip()
                interaction_type = str(row.get("type", "")).strip()

                if not mirna_name or not cerna_name:
                    continue

                if interaction_type not in PAIRED_COHORT_VALID_CERNA_EDGE_TYPES:
                    continue

                source_key = self.get_or_create_rna_file_node(
                    nodes=nodes,
                    node_name_index=node_name_index,
                    name=mirna_name,
                    rna_type="miRNA",
                    species=row.get("species"),
                    source="paired_cohort_cerna_axis",
                )

                target_rna_type = (
                    PAIRED_COHORT_CERNA_TYPE_TO_TARGET_RNA_TYPE[
                        interaction_type
                    ]
                )

                target_key = self.get_or_create_rna_file_node(
                    nodes=nodes,
                    node_name_index=node_name_index,
                    name=cerna_name,
                    rna_type=target_rna_type,
                    species=row.get("species"),
                    source="paired_cohort_cerna_axis",
                )

                pair_key = edge_pair_key(source_key, target_key)

                if pair_key in edge_pair_index:
                    existing_edge_id = edge_pair_index[pair_key]

                    append_database_to_edge(
                        edge=edges[existing_edge_id],
                        database=row.get("database"),
                    )

                    self.append_cerna_axis_info_to_edge(
                        edge=edges[existing_edge_id],
                        row=row,
                    )

                    continue

                edge_id = paired_cohort_cerna_edge_id(
                    source_name=mirna_name,
                    target_name=cerna_name,
                    edge_type=interaction_type,
                )

                edges[edge_id] = {
                    "id": edge_id,
                    "source": source_key,
                    "target": target_key,

                    # 对齐 CustomListQueryTaskNetworkView 的 edge contract。
                    "source_db_id": None,
                    "target_db_id": None,

                    "source_name": mirna_name,
                    "target_name": cerna_name,
                    "source_type": "miRNA",
                    "target_type": target_rna_type,
                    "species": str(row.get("species", "") or "").strip(),

                    # 保持和 Custom List 主网络边一致：
                    # type 表示相互作用类型。
                    "type": interaction_type,

                    # 关键：主网络边统一为 rna_interaction，
                    # 前端无需区分 cerna_axis / rna_interaction 两套样式。
                    "edge_type": "rna_interaction",

                    # Paired Cohort 专属来源。
                    "edge_source": "paired_cohort_cerna_axis",

                    "databases": parse_database_list(row.get("database")),
                    "is_immune_related": False,
                    "immune_annotation": None,

                    # Paired Cohort 特有字段。
                    "paired_cohort_annotation": {
                        "inference": str(
                            row.get("inference", "") or ""
                        ).strip(),
                        "miRNA_log2FC": to_float_or_none(
                            row.get("miRNA_log2FC")
                        ),
                        "ceRNA_log2FC": to_float_or_none(
                            row.get("ceRNA_log2FC")
                        ),
                    },

                    # 为兼容已有前端读取方式，也可以平铺保留。
                    "inference": str(row.get("inference", "") or "").strip(),
                    "miRNA_log2FC": to_float_or_none(
                        row.get("miRNA_log2FC")
                    ),
                    "ceRNA_log2FC": to_float_or_none(
                        row.get("ceRNA_log2FC")
                    ),
                }

                edge_pair_index[pair_key] = edge_id
                edge_count += 1

        return {
            "loaded": True,
            "edge_count": edge_count,
            "file": str(ceRNA_axis_file),
        }

    def add_immune_annotation_edges(
            self,
            immune_axis_file: Path,
            nodes: dict,
            edges: dict,
            edge_pair_index: dict,
            node_name_index: dict,
    ) -> dict:
        if not immune_axis_file.exists() or not immune_axis_file.is_file():
            return {
                "loaded": False,
                "edge_count": 0,
                "merged_edge_count": 0,
                "new_edge_count": 0,
                "file": str(immune_axis_file),
            }

        total_immune_records = 0
        merged_edge_count = 0
        new_edge_count = 0

        with immune_axis_file.open("r", encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f)

            if not reader.fieldnames:
                return {
                    "loaded": False,
                    "edge_count": 0,
                    "merged_edge_count": 0,
                    "new_edge_count": 0,
                    "file": str(immune_axis_file),
                }

            missing_columns = (
                    PAIRED_COHORT_IMMUNE_AXIS_REQUIRED_COLUMNS
                    - set(reader.fieldnames)
            )

            if missing_columns:
                raise ValueError(
                    "Missing required immune axis columns: "
                    f"{sorted(missing_columns)}"
                )

            for row in reader:
                mirna_name = str(row.get("miRNA", "")).strip()
                immune_gene = str(
                    row.get("Immune checkpointGene", "")
                ).strip()

                if not mirna_name or not immune_gene:
                    continue

                source_key = self.get_or_create_rna_file_node(
                    nodes=nodes,
                    node_name_index=node_name_index,
                    name=mirna_name,
                    rna_type="miRNA",
                    species="Homo sapiens",
                    source="immune_annotation",
                )

                target_key = self.get_or_create_immune_target_node(
                    nodes=nodes,
                    node_name_index=node_name_index,
                    immune_gene=immune_gene,
                )

                mark_node_as_immune_related(
                    nodes=nodes,
                    node_key=source_key,
                    immune_node_type="immune_source_miRNA",
                    pathway=row.get("Immune checkpointPathway"),
                )

                mark_node_as_immune_related(
                    nodes=nodes,
                    node_key=target_key,
                    immune_node_type="immune_checkpoint_gene",
                    pathway=row.get("Immune checkpointPathway"),
                )

                evidence_item = format_immune_evidence_item(row)

                pair_key = edge_pair_key(source_key, target_key)

                if pair_key in edge_pair_index:
                    existing_edge_id = edge_pair_index[pair_key]

                    append_immune_annotation_to_edge(
                        edge=edges[existing_edge_id],
                        evidence_item=evidence_item,
                    )

                    merged_edge_count += 1
                else:
                    edge_id = immune_edge_id(
                        mirna_name=mirna_name,
                        immune_gene=immune_gene,
                    )

                    if edge_id not in edges:
                        edges[edge_id] = {
                            "id": edge_id,
                            "source": source_key,
                            "target": target_key,

                            # 对齐 CustomListQueryTaskNetworkView 的 edge contract。
                            "source_db_id": None,
                            "target_db_id": None,

                            "source_name": mirna_name,
                            "target_name": immune_gene,
                            "source_type": "miRNA",
                            "target_type": "immune_checkpoint_gene",
                            "species": "Homo sapiens",
                            "type": "immune_annotation",
                            "edge_type": "immune_annotation",
                            "edge_source": "paired_cohort_immune_axis",
                            "databases": ["immune_annotation"],
                            "is_immune_related": True,
                            "immune_annotation": {
                                "evidence_count": 0,
                                "evidence_items": [],
                            },

                            # 为保持 Paired Cohort edge 字段稳定。
                            "paired_cohort_annotation": None,
                            "inference": "",
                            "miRNA_log2FC": None,
                            "ceRNA_log2FC": None,
                        }

                        edge_pair_index[pair_key] = edge_id
                        new_edge_count += 1

                    append_immune_annotation_to_edge(
                        edge=edges[edge_id],
                        evidence_item=evidence_item,
                    )

                total_immune_records += 1

        return {
            "loaded": True,
            "edge_count": total_immune_records,
            "merged_edge_count": merged_edge_count,
            "new_edge_count": new_edge_count,
            "file": str(immune_axis_file),
        }

    def get_or_create_rna_file_node(
            self,
            nodes: dict,
            node_name_index: dict,
            name: str,
            rna_type: str,
            species=None,
            source: str = "paired_cohort_cerna_axis",
    ) -> str:
        node_key = paired_cohort_rna_file_node_key(
            name=name,
            rna_type=rna_type,
        )

        if node_key not in nodes:
            nodes[node_key] = {
                "id": node_key,
                "db_id": None,
                "name": name,
                "type": rna_type,
                "species": str(species or "").strip() or None,
                "source": source,
                "matched_in_database": False,
                "is_immune_related": False,
                "immune_node_type": None,
                "immune_pathways": [],
            }

        normalized_name = normalize_rna_name(name)

        if normalized_name not in node_name_index:
            node_name_index[normalized_name] = []

        if node_key not in node_name_index[normalized_name]:
            node_name_index[normalized_name].append(node_key)

        return node_key

    def get_or_create_immune_target_node(
            self,
            nodes: dict,
            node_name_index: dict,
            immune_gene: str,
    ) -> str:
        normalized_name = normalize_rna_name(immune_gene)

        existing_node_keys = node_name_index.get(normalized_name, [])

        if existing_node_keys:
            return self.choose_preferred_existing_node(
                nodes=nodes,
                node_keys=existing_node_keys,
            )

        node_key = immune_gene_node_key(immune_gene)

        if node_key not in nodes:
            nodes[node_key] = {
                "id": node_key,
                "db_id": None,
                "name": immune_gene,
                "type": "immune_checkpoint_gene",
                "species": "Homo sapiens",
                "source": "immune_annotation",
                "matched_in_database": False,
                "is_immune_related": True,
                "immune_node_type": "immune_checkpoint_gene",
                "immune_pathways": [],
            }

        if normalized_name not in node_name_index:
            node_name_index[normalized_name] = []

        if node_key not in node_name_index[normalized_name]:
            node_name_index[normalized_name].append(node_key)

        return node_key

    @staticmethod
    def choose_preferred_existing_node(
            nodes: dict,
            node_keys: list[str],
    ) -> str:
        preferred_types = [
            "mRNA",
            "lncRNA",
            "circRNA",
            "miRNA",
            "immune_checkpoint_gene",
        ]

        for preferred_type in preferred_types:
            for node_key in node_keys:
                node = nodes.get(node_key)

                if node and node.get("type") == preferred_type:
                    return node_key

        return node_keys[0]

    @staticmethod
    def append_cerna_axis_info_to_edge(edge: dict, row: dict) -> None:
        """
        当 ceRNA_axis.csv 中同一 source-target 重复出现时，
        保留主 edge，同时补充 database 和 Paired Cohort 注释信息。

        当前前端主要依赖单条 edge 渲染，因此不把重复记录拆成多条边。
        """

        if "paired_cohort_annotation" not in edge:
            edge["paired_cohort_annotation"] = {
                "inference": str(row.get("inference", "") or "").strip(),
                "miRNA_log2FC": to_float_or_none(row.get("miRNA_log2FC")),
                "ceRNA_log2FC": to_float_or_none(row.get("ceRNA_log2FC")),
            }

        if not edge.get("inference"):
            edge["inference"] = str(row.get("inference", "") or "").strip()

        if edge.get("miRNA_log2FC") is None:
            edge["miRNA_log2FC"] = to_float_or_none(
                row.get("miRNA_log2FC")
            )

        if edge.get("ceRNA_log2FC") is None:
            edge["ceRNA_log2FC"] = to_float_or_none(
                row.get("ceRNA_log2FC")
            )


class WorkflowTaskResultDownloadView(APIView):
    """
    Download workflow task result archive.

    Query params:
        taskUUID=<uuid>

    Rules:
        1. Only completed Success tasks can be downloaded.
        2. The result archive is generated on demand.
        3. The generated zip archive is cached under task output directory.
        4. Subsequent downloads reuse the cached archive.
    """

    def get(self, request):
        try:
            task_uuid = str(
                request.query_params.get("taskUUID", "")
            ).strip()

            if not task_uuid:
                return Response(
                    {
                        "success": False,
                        "msg": "Missing required parameter: taskUUID.",
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                download_result = prepare_task_result_download(
                    task_uuid=task_uuid,
                )

            except InvalidTaskUUIDError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            except TaskNotReadyForDownloadError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            except TaskResultFileNotFoundError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            except TaskResultArchiveError as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_500_INTERNAL_SERVER_ERROR,
                )

            except (
                    TaskDownloadError,
                    TaskDownloadConfigNotFoundError,
                    TaskNotFoundError,
                    MultipleTaskMatchedError,
            ) as e:
                return Response(
                    {
                        "success": False,
                        "msg": str(e),
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            archive_path = download_result.archive_path

            if not archive_path.exists() or not archive_path.is_file():
                return Response(
                    {
                        "success": False,
                        "msg": "Task result archive not found.",
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            return FileResponse(
                open(archive_path, "rb"),
                as_attachment=True,
                filename=download_result.archive_name,
                content_type="application/zip",
            )

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {
                    "success": False,
                    "msg": f"Server error: {str(e)}",
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )


class PairedCohortAxisFinalDataView(APIView):
    """
    Return all rows from paired cohort ceRNA axis final result file.

    Query params:
        taskUUID: PairedCohortTask UUID

    Input filename:
        {task_name}_ceRNA_axis_final.csv
    """

    NUMERIC_COLUMNS = {
        "mRNA_log2FC",
        "miRNA_log2FC",
        "lncRNA_log2FC",
        "circRNA_log2FC",
    }

    def get(self, request):
        try:
            task_uuid = str(request.query_params.get("taskUUID", "")).strip()

            if not task_uuid:
                return Response(
                    {"detail": "Missing query parameter: taskUUID."},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                uuid_lib.UUID(task_uuid)
            except ValueError:
                return Response(
                    {"detail": "Invalid Task UUID format."},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                task = PairedCohortTask.objects.get(uuid=task_uuid)
            except PairedCohortTask.DoesNotExist:
                return Response(
                    {
                        "detail": f"PairedCohortTask not found: {task_uuid}."
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            if task.status != PairedCohortTask.Status.Success:
                return Response(
                    {
                        "detail": (
                            "Paired cohort task is not completed successfully. "
                            f"Current status: {task.get_status_display()}."
                        )
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                axis_file, df = read_paired_cohort_axis_final_file(task)
            except PairedCohortTaskPathError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_400_BAD_REQUEST,
                )
            except FileNotFoundError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_404_NOT_FOUND,
                )
            except PairedCohortTaskInputError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            df = normalize_paired_cohort_axis_final_dataframe(df)

            raw_count = int(df.shape[0])

            results = self.serialize_dataframe(df)

            return Response(
                {
                    "uuid": str(task.uuid),
                    "task_name": task.task_name,
                    "axis_final_file": axis_file.name,
                    "count": raw_count,
                    "columns": PAIRED_COHORT_AXIS_FINAL_COLUMNS,
                    "results": results,
                },
                status=status.HTTP_200_OK,
            )

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {"detail": f"Server error: {str(e)}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

    def serialize_dataframe(self, df: pd.DataFrame) -> list[dict]:
        clean_df = df.replace([np.inf, -np.inf], np.nan)

        records = []

        for _, row in clean_df.iterrows():
            item = {}

            for col in PAIRED_COHORT_AXIS_FINAL_COLUMNS:
                value = row.get(col)

                if col in self.NUMERIC_COLUMNS:
                    item[col] = self.safe_float_or_none(value)
                else:
                    item[col] = self.safe_str_or_empty(value)

            records.append(item)

        return records

    @staticmethod
    def safe_float_or_none(value):
        if value is None or pd.isna(value):
            return None

        try:
            result = float(value)
        except (TypeError, ValueError):
            return None

        if pd.isna(result) or not np.isfinite(result):
            return None

        return result

    @staticmethod
    def safe_str_or_empty(value):
        if value is None or pd.isna(value):
            return ""

        return str(value)


class PairedCohortSurvivalKMDataView(APIView):
    """
    Return Kaplan-Meier survival curve data for paired cohort task.

    Query params:
        taskUUID: PairedCohortTask UUID

    Input filename:
        {task_name}_survival_analysis.csv
    """

    def get(self, request):
        try:
            task_uuid = str(request.query_params.get("taskUUID", "")).strip()

            if not task_uuid:
                return Response(
                    {"detail": "Missing query parameter: taskUUID."},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                uuid_lib.UUID(task_uuid)
            except ValueError:
                return Response(
                    {"detail": "Invalid Task UUID format."},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                task = PairedCohortTask.objects.get(uuid=task_uuid)
            except PairedCohortTask.DoesNotExist:
                return Response(
                    {
                        "detail": f"PairedCohortTask not found: {task_uuid}."
                    },
                    status=status.HTTP_404_NOT_FOUND,
                )

            if task.status != PairedCohortTask.Status.Success:
                return Response(
                    {
                        "detail": (
                            "Paired cohort task is not completed successfully. "
                            f"Current status: {task.get_status_display()}."
                        )
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            try:
                result = build_paired_cohort_survival_km_data(task)
            except PairedCohortTaskPathError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_400_BAD_REQUEST,
                )
            except FileNotFoundError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_404_NOT_FOUND,
                )
            except PairedCohortTaskInputError as e:
                return Response(
                    {"detail": str(e)},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            return Response(result, status=status.HTTP_200_OK)

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {"detail": f"Server error: {str(e)}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )
