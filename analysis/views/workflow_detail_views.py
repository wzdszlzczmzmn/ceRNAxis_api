import csv
import traceback
import uuid as uuid_lib
from pathlib import Path

from django.db.models import Q

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from database.models import RNANode, RNAInteraction
from analysis.models import CustomListQueryTask

from analysis.utils.custom_list_query_task_utils import (
    get_task_output_dir,
    get_task_out_prefix,
)


class WorkflowRNAInteractionNetworkView(APIView):
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
