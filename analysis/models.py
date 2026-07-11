import os
import uuid

from django.conf import settings
from django.db import models


def default_rnas():
    return {
        "miRNA": [],
        "mRNA": [],
        "mRNA_up": [],
        "mRNA_down": [],
        "lncRNA": [],
        "circRNA": [],
    }


class CustomListQueryTask(models.Model):
    class Status(models.TextChoices):
        Success = "S", "Success"
        Pending = "P", "Pending"
        Failed = "F", "Failed"
        Running = "R", "Running"

    uuid = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
        editable=False,
    )

    user = models.CharField(
        max_length=300,
        blank=True,
        null=True,
    )

    task_name = models.CharField(
        max_length=200,
    )

    status = models.CharField(
        max_length=1,
        choices=Status.choices,
        default=Status.Pending,
    )

    create_time = models.DateTimeField(
        auto_now_add=True,
    )

    finish_time = models.DateTimeField(
        null=True,
        blank=True,
    )

    # 兼容旧任务：旧版本用于保存 ImmiRImmiR_ACC 这类 immune annotation file value。
    # 新 Module 1 不再依赖该字段。
    map_info = models.CharField(
        max_length=200,
        blank=True,
        default="",
        help_text=(
            "Deprecated for CustomListQueryTask. "
            "Previously selected immune annotation file value, e.g. ImmiRImmiR_ACC."
        ),
    )

    # 新脚本 run_module1.sh 使用的 cancer_type，例如 ACC、BRCA、LUAD。
    cancer_type = models.CharField(
        max_length=50,
        blank=True,
        default="",
        db_index=True,
        help_text="Cancer type used by run_module1.sh, e.g. ACC.",
    )

    has_mrna_direction = models.BooleanField(
        default=False,
        db_index=True,
        help_text=(
            "Whether mRNA input is directional. "
            "When true, mRNA_up and mRNA_down are used "
            "instead of the standard mRNA list."
        ),
    )

    rnas = models.JSONField(
        default=default_rnas,
        help_text=(
            "Input RNA lists grouped by miRNA, mRNA, "
            "mRNA_up, mRNA_down, lncRNA and circRNA."
        ),
    )

    class Meta:
        db_table = "custom_list_query_task"
        ordering = ["-create_time"]
        indexes = [
            models.Index(fields=["status"]),
            models.Index(fields=["cancer_type"]),
            models.Index(fields=["create_time"]),
        ]

    def __str__(self):
        return f"{self.task_name} ({self.uuid})"

    def get_workspace_dir_absolute_path(self):
        return os.path.join(
            settings.WORKSPACE_HOME,
            str(self.uuid),
        )

    def get_input_dir_absolute_path(self):
        return os.path.join(
            self.get_workspace_dir_absolute_path(),
            "input",
        )

    def get_output_dir_absolute_path(self):
        return os.path.join(
            self.get_workspace_dir_absolute_path(),
            "output",
        )

    def get_rna_count(self, rna_type: str) -> int:
        if not isinstance(self.rnas, dict):
            return 0

        value = self.rnas.get(rna_type, [])

        if not isinstance(value, list):
            return 0

        return len(value)

    @property
    def miRNA_count(self):
        return self.get_rna_count("miRNA")

    @property
    def mRNA_count(self):
        return self.get_rna_count("mRNA")

    @property
    def mRNA_up_count(self):
        return self.get_rna_count("mRNA_up")

    @property
    def mRNA_down_count(self):
        return self.get_rna_count("mRNA_down")

    @property
    def lncRNA_count(self):
        return self.get_rna_count("lncRNA")

    @property
    def circRNA_count(self):
        return self.get_rna_count("circRNA")

    @property
    def total_rna_count(self):
        return (
                self.miRNA_count
                + self.mRNA_count
                + self.mRNA_up_count
                + self.mRNA_down_count
                + self.lncRNA_count
                + self.circRNA_count
        )


class PairedCohortTask(models.Model):
    class Status(models.IntegerChoices):
        Pending = 0, "Pending"
        Running = 1, "Running"
        Success = 2, "Success"
        Failed = 3, "Failed"

    uuid = models.UUIDField(
        default=uuid.uuid4,
        editable=False,
        unique=True,
        db_index=True,
    )

    user = models.CharField(
        max_length=128,
        blank=True,
        default="",
    )

    task_name = models.CharField(
        max_length=255,
    )

    status = models.IntegerField(
        choices=Status.choices,
        default=Status.Pending,
    )

    map_info = models.CharField(
        max_length=255,
    )

    deg_method = models.CharField(
        max_length=32,
    )

    cancer_type = models.CharField(
        max_length=50,
        blank=True,
        default="",
    )

    use_padj = models.BooleanField(default=True)

    mrna_file = models.CharField(
        max_length=255,
        blank=True,
        default="",
    )

    mirna_file = models.CharField(
        max_length=255,
        blank=True,
        default="",
    )

    lncrna_file = models.CharField(
        max_length=255,
        blank=True,
        default="",
    )

    circrna_file = models.CharField(
        max_length=255,
        blank=True,
        default="",
    )

    meta_file = models.CharField(
        max_length=255,
        blank=True,
        default="",
    )

    logfc_cutoff_mrna = models.FloatField(default=1)
    padj_cutoff_mrna = models.FloatField(default=0.05)

    logfc_cutoff_mirna = models.FloatField(default=1)
    padj_cutoff_mirna = models.FloatField(default=0.05)

    logfc_cutoff_lncrna = models.FloatField(default=1)
    padj_cutoff_lncrna = models.FloatField(default=0.05)

    logfc_cutoff_circrna = models.FloatField(default=1.0)
    padj_cutoff_circrna = models.FloatField(default=0.1)

    create_time = models.DateTimeField(
        auto_now_add=True,
    )

    finish_time = models.DateTimeField(
        null=True,
        blank=True,
    )

    class Meta:
        db_table = "paired_cohort_task"
        ordering = ["-create_time"]

    def get_workspace_dir_absolute_path(self):
        return os.path.join(
            settings.WORKSPACE_HOME,
            str(self.uuid),
        )

    def get_input_dir_absolute_path(self):
        return os.path.join(
            self.get_workspace_dir_absolute_path(),
            "input",
        )

    def get_output_dir_absolute_path(self):
        return os.path.join(
            self.get_workspace_dir_absolute_path(),
            "output",
        )


class HybridReferenceTask(models.Model):
    class Status(models.IntegerChoices):
        Pending = 0, "Pending"
        Running = 1, "Running"
        Success = 2, "Success"
        Failed = 3, "Failed"

    class DEGMethod(models.TextChoices):
        limma = "limma", "limma"
        deseq2 = "deseq2", "deseq2"

    class LncRNAType(models.TextChoices):
        log2count = "log2count", "log2count"
        log2fpkm = "log2fpkm", "log2fpkm"
        log2fpkmuq = "log2fpkmuq", "log2fpkmuq"
        log2tpm = "log2tpm", "log2tpm"

    uuid = models.UUIDField(
        default=uuid.uuid4,
        editable=False,
        unique=True,
        db_index=True,
    )

    user = models.CharField(
        max_length=128,
        blank=True,
        default="",
    )

    task_name = models.CharField(
        max_length=255,
    )

    status = models.IntegerField(
        choices=Status.choices,
        default=Status.Pending,
    )

    tcga_type = models.CharField(
        max_length=64,
        db_index=True,
    )

    lncrna_type = models.CharField(
        max_length=32,
        choices=LncRNAType.choices,
        default=LncRNAType.log2tpm,
    )

    use_padj = models.BooleanField(
        default=True,
    )

    map_info = models.CharField(
        max_length=255,
    )

    deg_method = models.CharField(
        max_length=32,
        choices=DEGMethod.choices,
        default=DEGMethod.limma,
    )

    mrna_file = models.CharField(
        max_length=255,
        blank=True,
        default="",
    )

    meta_file = models.CharField(
        max_length=255,
        blank=True,
        default="",
    )

    logfc_cutoff_mrna = models.FloatField(default=1.0)
    padj_cutoff_mrna = models.FloatField(default=0.05)

    create_time = models.DateTimeField(
        auto_now_add=True,
    )

    finish_time = models.DateTimeField(
        null=True,
        blank=True,
    )

    class Meta:
        db_table = "hybrid_reference_task"
        ordering = ["-create_time"]
        indexes = [
            models.Index(fields=["status"]),
            models.Index(fields=["tcga_type"]),
            models.Index(fields=["create_time"]),
        ]

    def __str__(self):
        return f"{self.task_name} ({self.uuid})"

    def get_workspace_dir_absolute_path(self):
        return os.path.join(
            settings.WORKSPACE_HOME,
            str(self.uuid),
        )

    def get_input_dir_absolute_path(self):
        return os.path.join(
            self.get_workspace_dir_absolute_path(),
            "input",
        )

    def get_output_dir_absolute_path(self):
        return os.path.join(
            self.get_workspace_dir_absolute_path(),
            "output",
        )
