import csv
from pathlib import Path

from analysis.utils.paired_cohort_network_utils import (
    paired_cohort_rna_file_node_key,
    PAIRED_COHORT_IMMUNE_AXIS_REQUIRED_COLUMNS,
    parse_database_list,
    paired_cohort_cerna_edge_id,
    append_database_to_edge,
    PAIRED_COHORT_CERNA_TYPE_TO_TARGET_RNA_TYPE,
    PAIRED_COHORT_VALID_CERNA_EDGE_TYPES,
    PAIRED_COHORT_CERNA_AXIS_REQUIRED_COLUMNS,
)
from analysis.utils.workflow_detail_utils.workflow_network_utils import (
    immune_gene_node_key,
    normalize_rna_name,
    append_immune_annotation_to_edge,
    immune_edge_id,
    edge_pair_key,
    format_immune_evidence_item,
    mark_node_as_immune_related,
    to_float_or_none,
)


class WorkflowOutputNetworkBuilderMixin:
    cerna_node_source = "workflow_cerna_axis"
    cerna_edge_source = "workflow_cerna_axis"
    immune_node_source = "immune_annotation"
    immune_edge_source = "workflow_immune_axis"

    def build_network_from_files(
        self,
        *,
        ceRNA_axis_file: Path,
        immune_axis_file: Path,
    ) -> dict:
        if not ceRNA_axis_file.exists() or not ceRNA_axis_file.is_file():
            return self.empty_network_response(
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

                "rna_interaction_edge_count": cerna_result["edge_count"],

                "cerna_axis_edge_count": cerna_result["edge_count"],
                "cerna_axis_loaded": cerna_result["loaded"],

                "immune_annotation_edge_count": immune_result["edge_count"],
                "immune_merged_edge_count": immune_result["merged_edge_count"],
                "immune_new_edge_count": immune_result["new_edge_count"],
                "immune_annotation_loaded": immune_result["loaded"],

                "cerna_axis_file": cerna_result["file"],
                "immune_annotation_file": immune_result["file"],

                "edge_scope": "dataset_annotation_output_files",
            },
            "nodes": list(nodes.values()),
            "edges": list(edges.values()),
            "ignored_nodes": [],
        }

    def empty_network_response(
        self,
        ceRNA_axis_file: Path,
        immune_axis_file: Path,
    ) -> dict:
        return {
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
                "edge_scope": "dataset_annotation_output_files",
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
                    source=self.cerna_node_source,
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
                    source=self.cerna_node_source,
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

                    "source_db_id": None,
                    "target_db_id": None,

                    "source_name": mirna_name,
                    "target_name": cerna_name,
                    "source_type": "miRNA",
                    "target_type": target_rna_type,
                    "species": str(row.get("species", "") or "").strip(),

                    "type": interaction_type,
                    "edge_type": "rna_interaction",
                    "edge_source": self.cerna_edge_source,

                    "databases": parse_database_list(row.get("database")),
                    "is_immune_related": False,
                    "immune_annotation": None,

                    "workflow_annotation": {
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
                    source=self.immune_node_source,
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

                            "source_db_id": None,
                            "target_db_id": None,

                            "source_name": mirna_name,
                            "target_name": immune_gene,
                            "source_type": "miRNA",
                            "target_type": "immune_checkpoint_gene",
                            "species": "Homo sapiens",
                            "type": "immune_annotation",
                            "edge_type": "immune_annotation",
                            "edge_source": self.immune_edge_source,
                            "databases": ["immune_annotation"],
                            "is_immune_related": True,
                            "immune_annotation": {
                                "evidence_count": 0,
                                "evidence_items": [],
                            },

                            "workflow_annotation": None,
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
        source: str = "workflow_cerna_axis",
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
                "source": self.immune_node_source,
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
        if "workflow_annotation" not in edge:
            edge["workflow_annotation"] = {
                "inference": str(row.get("inference", "") or "").strip(),
                "miRNA_log2FC": to_float_or_none(row.get("miRNA_log2FC")),
                "ceRNA_log2FC": to_float_or_none(row.get("ceRNA_log2FC")),
            }

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
