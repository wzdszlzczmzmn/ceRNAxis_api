import os
import uuid

from django.conf import settings
from django.db import models


def default_rnas():
    return {
        "miRNA": [],
        "mRNA": [],
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

    map_info = models.CharField(
        max_length=200,
        help_text="Selected immune annotation file value, e.g. ImmiRImmiR_ACC.",
    )

    rnas = models.JSONField(
        default=default_rnas,
        help_text="Input RNA lists grouped by miRNA, mRNA, lncRNA and circRNA.",
    )

    class Meta:
        db_table = "custom_list_query_task"
        ordering = ["-create_time"]

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
