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
