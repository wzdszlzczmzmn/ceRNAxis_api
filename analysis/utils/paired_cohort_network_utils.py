from pathlib import Path

from analysis.utils.paired_cohort_task_utils import validate_safe_name, get_paired_cohort_task_output_dir, \
    PairedCohortTaskPathError

PAIRED_COHORT_VALID_CERNA_EDGE_TYPES = {
    "miRNA-mRNA",
    "miRNA-lncRNA",
    "miRNA-circRNA",
}

PAIRED_COHORT_CERNA_TYPE_TO_TARGET_RNA_TYPE = {
    "miRNA-mRNA": "mRNA",
    "miRNA-lncRNA": "lncRNA",
    "miRNA-circRNA": "circRNA",
}

PAIRED_COHORT_CERNA_AXIS_REQUIRED_COLUMNS = {
    "miRNA",
    "ceRNA",
    "species",
    "database",
    "type",
    "miRNA_log2FC",
    "ceRNA_log2FC",
    "inference",
}

PAIRED_COHORT_IMMUNE_AXIS_REQUIRED_COLUMNS = {
    "miRNA",
    "Immune checkpointGene",
}


def get_paired_cohort_cerna_axis_file_path(task) -> Path:
    task_name = str(task.task_name).strip()

    validate_safe_name(task_name, "task_name")

    output_dir = get_paired_cohort_task_output_dir(task)
    file_path = (output_dir / f"{task_name}_ceRNA_axis.csv").resolve()

    if not str(file_path).startswith(str(output_dir)):
        raise PairedCohortTaskPathError(
            "Invalid paired cohort ceRNA axis file path."
        )

    return file_path


def get_paired_cohort_immune_axis_file_path(task) -> Path:
    task_name = str(task.task_name).strip()

    validate_safe_name(task_name, "task_name")

    output_dir = get_paired_cohort_task_output_dir(task)
    file_path = (output_dir / f"{task_name}_map_immune_axis.csv").resolve()

    if not str(file_path).startswith(str(output_dir)):
        raise PairedCohortTaskPathError(
            "Invalid paired cohort immune axis file path."
        )

    return file_path


def paired_cohort_rna_file_node_key(name: str, rna_type: str) -> str:
    return f"rna:{rna_type}:{name}"


def paired_cohort_cerna_edge_id(
    source_name: str,
    target_name: str,
    edge_type: str,
) -> str:
    return f"cerna_axis:{edge_type}:{source_name}:{target_name}"


def parse_database_list(database) -> list[str]:
    value = str(database or "").strip()

    if not value:
        return []

    parts = []

    for sep in [";", "|", ","]:
        if sep in value:
            parts = [item.strip() for item in value.split(sep)]
            break

    if not parts:
        parts = [value]

    return [
        item
        for item in parts
        if item
    ]


def append_database_to_edge(edge: dict, database) -> None:
    databases = edge.setdefault("databases", [])

    for db in parse_database_list(database):
        if db not in databases:
            databases.append(db)
