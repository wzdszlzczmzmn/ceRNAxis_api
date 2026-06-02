import json
import random
from pathlib import Path

from django.core.management.base import BaseCommand
from django.db.models import Q

from database.models import RNANode, RNAInteraction


class Command(BaseCommand):
    help = "Export a structured JSON network from sampled miRNA/mRNA and predefined lncRNA/circRNA nodes."

    LNC_RNAS = [
        "MALAT1",
        "HOTAIR",
        "XIST",
        "NEAT1",
        "H19",
    ]

    CIRC_RNAS = [
        "CDR1as",
        "circHIPK3",
        "circSMARCA5",
        "circITCH",
        "circFOXO3",
    ]

    def add_arguments(self, parser):
        parser.add_argument(
            "--output",
            type=str,
            default=None,
            help="Optional output JSON file path.",
        )

        parser.add_argument(
            "--sample-size",
            type=int,
            default=5,
            help="Random sample size for miRNA and mRNA.",
        )

    def handle(self, *args, **options):
        sample_size = options["sample_size"]
        output_path = options["output"]

        sampled_mirnas = self.sample_nodes("miRNA", sample_size)
        sampled_mrnas = self.sample_nodes("mRNA", sample_size)

        seed_items = []

        for node in sampled_mirnas:
            seed_items.append({
                "name": node.name,
                "type": "miRNA",
            })

        for node in sampled_mrnas:
            seed_items.append({
                "name": node.name,
                "type": "mRNA",
            })

        for name in self.LNC_RNAS:
            seed_items.append({
                "name": name,
                "type": "lncRNA",
            })

        for name in self.CIRC_RNAS:
            seed_items.append({
                "name": name,
                "type": "circRNA",
            })

        result = self.build_network(seed_items)

        json_text = json.dumps(
            result,
            ensure_ascii=False,
            indent=2,
        )

        if output_path:
            path = Path(output_path)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json_text, encoding="utf-8")
            self.stdout.write(
                self.style.SUCCESS(f"Exported JSON to {path}")
            )
        else:
            self.stdout.write(json_text)

    def sample_nodes(self, rna_type, sample_size):
        queryset = RNANode.objects.filter(rna_type=rna_type)

        count = queryset.count()

        if count == 0:
            return []

        sample_size = min(sample_size, count)
        offsets = random.sample(range(count), sample_size)

        nodes = []

        for offset in offsets:
            node = queryset.order_by("id")[offset]
            nodes.append(node)

        return nodes

    def build_network(self, seed_items):
        seed_key_set = {
            (item["name"], item["type"])
            for item in seed_items
        }

        node_query = Q()

        for name, rna_type in seed_key_set:
            node_query |= Q(name=name, rna_type=rna_type)

        existing_nodes = []

        if node_query:
            existing_nodes = list(
                RNANode.objects.filter(node_query)
            )

        existing_node_map = {
            (node.name, node.rna_type): node
            for node in existing_nodes
        }

        seed_db_ids = [
            node.id
            for node in existing_nodes
        ]

        nodes = {}
        edges = {}

        # 先加入所有输入节点：即使数据库不存在，也保留
        for item in seed_items:
            name = item["name"]
            rna_type = item["type"]

            db_node = existing_node_map.get((name, rna_type))
            key = self.node_key(db_node.id if db_node else None, name, rna_type)

            nodes[key] = {
                "id": db_node.id if db_node else None,
                "name": name,
                "type": rna_type,
                "species": db_node.species if db_node else None,
                "in_database": db_node is not None,
                "is_seed": True,
            }

        # 关键变化：
        # 只保留给定节点列表内部的互作关系
        if seed_db_ids:
            interactions = (
                RNAInteraction.objects
                .filter(
                    source_id__in=seed_db_ids,
                    target_id__in=seed_db_ids,
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

            source_key = self.node_key(source.id, source.name, source.rna_type)
            target_key = self.node_key(target.id, target.name, target.rna_type)

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
                "source_name": source.name,
                "target_name": target.name,
                "source_type": source.rna_type,
                "target_type": target.rna_type,
                "species": interaction.species,
                "type": interaction.interaction_type,
                "databases": [
                    db.name for db in interaction.databases.all()
                ],
            }

        return {
            "meta": {
                "seed_node_count": len(seed_items),
                "database_seed_node_count": len(existing_nodes),
                "total_node_count": len(nodes),
                "total_edge_count": len(edges),
                "edge_scope": "seed_nodes_only",
            },
            "seed_nodes": seed_items,
            "nodes": list(nodes.values()),
            "edges": list(edges.values()),
        }

    def node_key(self, node_id, name, rna_type):
        if node_id is not None:
            return f"node:{node_id}"

        return f"external:{rna_type}:{name}"
