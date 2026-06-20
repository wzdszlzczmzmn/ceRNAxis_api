import math
import uuid as uuid_lib


def is_valid_uuid(value: str) -> bool:
    try:
        uuid_lib.UUID(str(value))
        return True
    except ValueError:
        return False


def normalize_rna_name(name: str) -> str:
    return str(name or "").strip().lower()


def to_float_or_none(value):
    if value is None:
        return None

    value = str(value).strip()

    if not value:
        return None

    try:
        result = float(value)
    except ValueError:
        return None

    if math.isnan(result) or math.isinf(result):
        return None

    return result


def format_immune_evidence_item(row: dict) -> dict:
    return {
        "correlation_coefficient": to_float_or_none(
            row.get("CorrelationCoefficient")
        ),
        "p_value": to_float_or_none(
            row.get("P-value")
        ),
        "q_value": to_float_or_none(
            row.get("Q-value")
        ),
        "cancer": str(row.get("Cancer", "") or "").strip(),
        "pathway": str(row.get("Immune checkpointPathway", "") or "").strip(),
        "evidence": str(row.get("Evidence", "") or "").strip(),
    }


def immune_evidence_key(evidence_item: dict):
    return (
        evidence_item.get("correlation_coefficient"),
        evidence_item.get("p_value"),
        evidence_item.get("q_value"),
        evidence_item.get("cancer"),
        evidence_item.get("pathway"),
        evidence_item.get("evidence"),
    )


def append_immune_annotation_to_edge(
    edge: dict,
    evidence_item: dict,
) -> None:
    edge["is_immune_related"] = True

    if edge.get("immune_annotation") is None:
        edge["immune_annotation"] = {
            "evidence_count": 0,
            "evidence_items": [],
        }

    evidence_items = edge["immune_annotation"]["evidence_items"]

    evidence_key = immune_evidence_key(evidence_item)

    existing_keys = {
        immune_evidence_key(item)
        for item in evidence_items
    }

    if evidence_key not in existing_keys:
        evidence_items.append(evidence_item)

    edge["immune_annotation"]["evidence_count"] = len(evidence_items)


def mark_node_as_immune_related(
    nodes: dict,
    node_key: str,
    immune_node_type: str,
    pathway=None,
) -> None:
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


def edge_pair_key(source_key: str, target_key: str) -> str:
    return f"{source_key}=>{target_key}"


def immune_gene_node_key(gene_name: str) -> str:
    return f"immune_checkpoint_gene:{gene_name}"


def immune_edge_id(mirna_name: str, immune_gene: str) -> str:
    return f"immune_annotation:{mirna_name}:{immune_gene}"
