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
