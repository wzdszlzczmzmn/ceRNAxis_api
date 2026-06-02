from django.db.models import Q
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from database.models import RNANode, RNAInteraction


class RNAInteractionNetworkView(APIView):
    """
    根据前端传入的 RNA 列表查询网络。

    输入示例：
    {
        "miRNA": [...],
        "mRNA": [...],
        "lncRNA": [...],
        "circRNA": [...]
    }

    查询语义：
    1. 只保留数据库中存在的节点
    2. 不存在的节点放入 ignored_nodes
    3. 只查询这些已存在节点之间的互作关系
    4. 不展开节点的所有外部互作
    """

    ALLOWED_RNA_TYPES = {
        "miRNA",
        "mRNA",
        "lncRNA",
        "circRNA",
    }

    def post(self, request, *args, **kwargs):
        input_items = self.parse_input(request.data)

        if not input_items:
            return Response(
                {
                    "meta": {
                        "input_node_count": 0,
                        "matched_node_count": 0,
                        "ignored_node_count": 0,
                        "edge_count": 0,
                    },
                    "nodes": [],
                    "edges": [],
                    "ignored_nodes": [],
                },
                status=status.HTTP_200_OK,
            )

        result = self.build_network(input_items)

        return Response(result, status=status.HTTP_200_OK)

    def parse_input(self, data):
        """
        将前端分组结构转成：
        [
            {"name": "MALAT1", "type": "lncRNA"},
            ...
        ]
        """
        items = []
        seen = set()

        for rna_type in self.ALLOWED_RNA_TYPES:
            names = data.get(rna_type, [])

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

                items.append({
                    "name": name,
                    "type": rna_type,
                })

        return items

    def build_network(self, input_items):
        input_key_set = {
            (item["name"], item["type"])
            for item in input_items
        }

        node_query = Q()

        for name, rna_type in input_key_set:
            node_query |= Q(
                name=name,
                rna_type=rna_type,
            )

        matched_nodes = list(
            RNANode.objects.filter(node_query)
        )

        matched_node_map = {
            (node.name, node.rna_type): node
            for node in matched_nodes
        }

        matched_node_ids = [
            node.id
            for node in matched_nodes
        ]

        ignored_nodes = []

        for item in input_items:
            key = (item["name"], item["type"])

            if key not in matched_node_map:
                ignored_nodes.append({
                    "name": item["name"],
                    "type": item["type"],
                    "reason": "not_found_in_database",
                })

        nodes = {}
        edges = {}
        used_node_keys = set()

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

        for interaction in interactions:
            source = interaction.source
            target = interaction.target

            source_key = self.node_key(source)
            target_key = self.node_key(target)

            used_node_keys.add(source_key)
            used_node_keys.add(target_key)

            if source_key not in nodes:
                nodes[source_key] = {
                    "id": source.id,
                    "name": source.name,
                    "type": source.rna_type,
                    "species": source.species,
                }

            if target_key not in nodes:
                nodes[target_key] = {
                    "id": target.id,
                    "name": target.name,
                    "type": target.rna_type,
                    "species": target.species,
                }

            edge_key = (
                interaction.source_id,
                interaction.target_id,
                interaction.species,
                interaction.interaction_type,
            )

            edges[edge_key] = {
                "id": interaction.id,
                "source": source_key,
                "target": target_key,
                "source_id": source.id,
                "target_id": target.id,
                "source_name": source.name,
                "target_name": target.name,
                "source_type": source.rna_type,
                "target_type": target.rna_type,
                "species": interaction.species,
                "type": interaction.interaction_type,
                "databases": [
                    db.name
                    for db in interaction.databases.all()
                ],
            }

        for node in matched_nodes:
            key = self.node_key(node)

            if key not in used_node_keys:
                ignored_nodes.append({
                    "name": node.name,
                    "type": node.rna_type,
                    "reason": "no_interaction_in_query_scope",
                })

        return {
            "meta": {
                "input_node_count": len(input_items),
                "matched_node_count": len(matched_nodes),
                "returned_node_count": len(nodes),
                "ignored_node_count": len(ignored_nodes),
                "edge_count": len(edges),
                "edge_scope": "matched_input_nodes_only",
            },
            "nodes": list(nodes.values()),
            "edges": list(edges.values()),
            "ignored_nodes": ignored_nodes,
        }

    def node_key(self, node):
        return f"node:{node.id}"
