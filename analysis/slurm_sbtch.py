import subprocess

from analysis.models import CustomListQueryTask, PairedCohortTask, HybridReferenceTask
from analysis.utils.custom_list_query_task_utils import get_task_output_dir
from analysis.utils.hybrid_reference_task_utils import validate_hybrid_reference_input_files, \
    get_hybrid_reference_task_output_dir, HYBRID_REFERENCE_VALID_TCGA_TYPES, HYBRID_REFERENCE_VALID_LNCRNA_TYPES, \
    HYBRID_REFERENCE_VALID_DEG_METHODS
from analysis.utils.immune_annotation_path_utils import validate_immune_annotation_file
from analysis.utils.paired_cohort_task_utils import validate_paired_cohort_input_files, \
    get_paired_cohort_task_output_dir
from analysis.utils.slurm_path_utils import validate_slurm_script, get_pipeline_stdout_file_path, \
    get_pipeline_stderr_file_path, get_slurm_job_name

CUSTOM_LIST_QUERY_SLURM_SCRIPT = "submit_custom_list_query_task.sh"
PAIRED_COHORT_SLURM_SCRIPT = "submit_paired_cohort_task.sh"
HYBRID_REFERENCE_SLURM_SCRIPT = "submit_hybrid_reference_task.sh"


def _join_rna_values(values) -> str:
    """Convert a normalized RNA list to a comma-separated string."""

    if not values:
        return ""

    if not isinstance(values, (list, tuple, set)):
        return ""

    cleaned_values = []

    for value in values:
        cleaned_value = str(value).strip()

        if cleaned_value:
            cleaned_values.append(cleaned_value)

    return ",".join(cleaned_values)


def sbatch_custom_list_query_task(task_uuid) -> dict:
    """
    Submit CustomListQueryTask to SLURM.

    Expected SLURM script arguments:
        1. uuid
        2. miRNA_str
        3. mRNA_str
        4. mRNA_str_up
        5. mRNA_str_down
        6. lncRNA_str
        7. circRNA_str
        8. outdir
        9. out_prefix
        10. cancer_type
        11. has_mRNA_direction
    """

    try:
        task = CustomListQueryTask.objects.get(uuid=task_uuid)
    except CustomListQueryTask.DoesNotExist:
        return {
            "success": False,
            "msg": f"CustomListQueryTask not found: {task_uuid}",
            "stdout": "",
            "stderr": "",
        }

    try:
        script_path = validate_slurm_script(
            CUSTOM_LIST_QUERY_SLURM_SCRIPT
        )

        output_dir = get_task_output_dir(task)
        output_dir.mkdir(parents=True, exist_ok=True)

        stdout_file = get_pipeline_stdout_file_path(output_dir)
        stderr_file = get_pipeline_stderr_file_path(output_dir)

        out_prefix = task.task_name

        cancer_type = str(getattr(task, "cancer_type", "") or "").strip()

        if not cancer_type:
            raise ValueError("Missing cancer_type for custom list query task.")

        rnas = task.rnas or {}

        miRNA_str = _join_rna_values(rnas.get("miRNA"))
        mRNA_str = _join_rna_values(rnas.get("mRNA"))
        lncRNA_str = _join_rna_values(rnas.get("lncRNA"))
        circRNA_str = _join_rna_values(rnas.get("circRNA"))

        # 当前阶段固定走非 directional mRNA 分支。
        mRNA_str_up = ""
        mRNA_str_down = ""
        has_mRNA_direction = "False"

        # 如果以后开放 directional 分支，可以改成：
        # has_mRNA_direction = "True" if task.has_mrna_direction else "False"

    except Exception as e:
        return {
            "success": False,
            "msg": str(e),
            "stdout": "",
            "stderr": "",
        }

    command = [
        "sbatch",
        f"--job-name={get_slurm_job_name(task.uuid)}",
        f"--output={stdout_file}",
        f"--error={stderr_file}",
        str(script_path),
        str(task.uuid),
        miRNA_str,
        mRNA_str,
        mRNA_str_up,
        mRNA_str_down,
        lncRNA_str,
        circRNA_str,
        str(output_dir),
        out_prefix,
        cancer_type,
        has_mRNA_direction,
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        return {
            "success": False,
            "msg": "Failed to submit SLURM job.",
            "stdout": result.stdout,
            "stderr": result.stderr,
            "command": command,
        }

    return {
        "success": True,
        "msg": "SLURM job submitted successfully.",
        "stdout": result.stdout,
        "stderr": result.stderr,
        "command": command,
    }


def sbatch_paired_cohort_task(task_uuid) -> dict:
    """
    Submit PairedCohortTask to SLURM.

    Expected SLURM script arguments:
        1. uuid
        2. dataset
        3. mrna_file
        4. mirna_file
        5. lncrna_file
        6. meta_file
        7. outdir
        8. logfc_cutoff_mrna
        9. padj_cutoff_mrna
        10. logfc_cutoff_mirna
        11. padj_cutoff_mirna
        12. logfc_cutoff_lncrna
        13. padj_cutoff_lncrna
        14. deg_method
        15. map_info_csv
    """

    try:
        task = PairedCohortTask.objects.get(uuid=task_uuid)
    except PairedCohortTask.DoesNotExist:
        return {
            "success": False,
            "msg": f"PairedCohortTask not found: {task_uuid}",
            "stdout": "",
            "stderr": "",
        }

    try:
        script_path = validate_slurm_script(
            PAIRED_COHORT_SLURM_SCRIPT
        )

        map_info_file = validate_immune_annotation_file(
            task.map_info
        )

        input_files = validate_paired_cohort_input_files(
            task
        )

        output_dir = get_paired_cohort_task_output_dir(task)
        output_dir.mkdir(parents=True, exist_ok=True)

        stdout_file = get_pipeline_stdout_file_path(output_dir)
        stderr_file = get_pipeline_stderr_file_path(output_dir)

        dataset = str(task.task_name).strip()

        if not dataset:
            raise ValueError("Missing task_name for paired cohort task.")

        if task.deg_method not in ["limma", "deseq2"]:
            raise ValueError(
                "Invalid deg_method. Allowed values are: limma, deseq2."
            )

        cancer_type = str(getattr(task, "cancer_type", "") or "").strip()
        use_padj = "TRUE" if getattr(task, "use_padj", True) else "FALSE"

    except Exception as e:
        return {
            "success": False,
            "msg": str(e),
            "stdout": "",
            "stderr": "",
        }

    command = [
        "sbatch",
        f"--job-name={get_slurm_job_name(task.uuid)}",
        f"--output={stdout_file}",
        f"--error={stderr_file}",
        str(script_path),
        str(task.uuid),
        dataset,
        str(input_files["mrna_file"]),
        str(input_files["mirna_file"]),
        str(input_files["lncrna_file"] or ""),
        str(input_files["circrna_file"] or ""),
        str(input_files["meta_file"]),
        str(output_dir),
        str(task.logfc_cutoff_mrna),
        str(task.padj_cutoff_mrna),
        str(task.logfc_cutoff_mirna),
        str(task.padj_cutoff_mirna),
        str(task.logfc_cutoff_lncrna),
        str(task.padj_cutoff_lncrna),
        str(task.logfc_cutoff_circrna),
        str(task.padj_cutoff_circrna),
        task.deg_method,
        str(map_info_file),
        cancer_type,
        use_padj,
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        return {
            "success": False,
            "msg": "Failed to submit SLURM job.",
            "stdout": result.stdout,
            "stderr": result.stderr,
            "command": command,
        }

    return {
        "success": True,
        "msg": "SLURM job submitted successfully.",
        "stdout": result.stdout,
        "stderr": result.stderr,
        "command": command,
    }


def sbatch_hybrid_reference_task(task_uuid) -> dict:
    """
    Submit HybridReferenceTask to SLURM.

    Expected SLURM script arguments:
        1. uuid
        2. dataset
        3. mrna_file
        4. meta_file
        5. tcga_type
        6. lncrna_type
        7. outdir
        8. logfc_cutoff_mrna
        9. padj_cutoff_mrna
        10. deg_method
        11. map_info_csv
        12. use_padj
    """

    try:
        task = HybridReferenceTask.objects.get(uuid=task_uuid)
    except HybridReferenceTask.DoesNotExist:
        return {
            "success": False,
            "msg": f"HybridReferenceTask not found: {task_uuid}",
            "stdout": "",
            "stderr": "",
        }

    try:
        script_path = validate_slurm_script(
            HYBRID_REFERENCE_SLURM_SCRIPT
        )

        map_info_file = validate_immune_annotation_file(
            task.map_info
        )

        input_files = validate_hybrid_reference_input_files(
            task
        )

        output_dir = get_hybrid_reference_task_output_dir(task)
        output_dir.mkdir(parents=True, exist_ok=True)

        stdout_file = get_pipeline_stdout_file_path(output_dir)
        stderr_file = get_pipeline_stderr_file_path(output_dir)

        dataset = str(task.task_name).strip()

        if not dataset:
            raise ValueError("Missing task_name for hybrid reference task.")

        if task.tcga_type not in HYBRID_REFERENCE_VALID_TCGA_TYPES:
            raise ValueError(
                "Invalid tcga_type. Allowed values are: "
                f"{', '.join(HYBRID_REFERENCE_VALID_TCGA_TYPES)}."
            )

        if task.lncrna_type not in HYBRID_REFERENCE_VALID_LNCRNA_TYPES:
            raise ValueError(
                "Invalid lncrna_type. Allowed values are: "
                f"{', '.join(HYBRID_REFERENCE_VALID_LNCRNA_TYPES)}."
            )

        if task.deg_method not in HYBRID_REFERENCE_VALID_DEG_METHODS:
            raise ValueError(
                "Invalid deg_method. Allowed values are: "
                f"{', '.join(HYBRID_REFERENCE_VALID_DEG_METHODS)}."
            )

        if not isinstance(task.use_padj, bool):
            raise ValueError(
                "Invalid use_padj. Allowed values are True or False."
            )

        use_padj = "TRUE" if task.use_padj else "FALSE"

    except Exception as e:
        return {
            "success": False,
            "msg": str(e),
            "stdout": "",
            "stderr": "",
        }

    command = [
        "sbatch",
        f"--job-name={get_slurm_job_name(task.uuid)}",
        f"--output={stdout_file}",
        f"--error={stderr_file}",
        str(script_path),
        str(task.uuid),
        dataset,
        str(input_files["mrna_file"]),
        str(input_files["meta_file"]),
        task.tcga_type,
        task.lncrna_type,
        str(output_dir),
        str(task.logfc_cutoff_mrna),
        str(task.padj_cutoff_mrna),
        task.deg_method,
        str(map_info_file),
        use_padj,
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        return {
            "success": False,
            "msg": "Failed to submit SLURM job.",
            "stdout": result.stdout,
            "stderr": result.stderr,
            "command": command,
        }

    return {
        "success": True,
        "msg": "SLURM job submitted successfully.",
        "stdout": result.stdout,
        "stderr": result.stderr,
        "command": command,
    }
