from django.utils import timezone


def format_datetime(value):
    if value is None:
        return None

    return timezone.localtime(value).strftime("%Y-%m-%d %H:%M:%S")


def normalize_task_rnas_for_response(task) -> dict:
    if not isinstance(task.rnas, dict):
        return {
            "miRNA": [],
            "mRNA": [],
            "lncRNA": [],
            "circRNA": [],
        }

    return {
        "miRNA": task.rnas.get("miRNA", []) if isinstance(task.rnas.get("miRNA", []), list) else [],
        "mRNA": task.rnas.get("mRNA", []) if isinstance(task.rnas.get("mRNA", []), list) else [],
        "lncRNA": task.rnas.get("lncRNA", []) if isinstance(task.rnas.get("lncRNA", []), list) else [],
        "circRNA": task.rnas.get("circRNA", []) if isinstance(task.rnas.get("circRNA", []), list) else [],
    }


def format_custom_list_query_task(task, position=0) -> dict:
    return {
        "uuid": str(task.uuid),
        "status": task.status,
        "position": position,
        "task_name": task.task_name,
        "map_info": task.map_info,
        "create_time": format_datetime(task.create_time),
        "finish_time": format_datetime(task.finish_time),
        "miRNA_count": task.miRNA_count,
        "mRNA_count": task.mRNA_count,
        "lncRNA_count": task.lncRNA_count,
        "circRNA_count": task.circRNA_count,
        "rnas": normalize_task_rnas_for_response(task),
    }


def format_paired_cohort_task(task, position=0) -> dict:
    return {
        "uuid": str(task.uuid),
        "status": task.status,
        "position": position,
        "task_name": task.task_name,
        "map_info": task.map_info,
        "deg_method": task.deg_method,
        "create_time": format_datetime(task.create_time),
        "finish_time": format_datetime(task.finish_time),
        "files": {
            "mrna_file": task.mrna_file,
            "mirna_file": task.mirna_file,
            "lncrna_file": task.lncrna_file,
            "meta_file": task.meta_file,
        },
        "cutoffs": {
            "mRNA": {
                "logfc_cutoff": task.logfc_cutoff_mrna,
                "padj_cutoff": task.padj_cutoff_mrna,
            },
            "miRNA": {
                "logfc_cutoff": task.logfc_cutoff_mirna,
                "padj_cutoff": task.padj_cutoff_mirna,
            },
            "lncRNA": {
                "logfc_cutoff": task.logfc_cutoff_lncrna,
                "padj_cutoff": task.padj_cutoff_lncrna,
            },
        },
    }
