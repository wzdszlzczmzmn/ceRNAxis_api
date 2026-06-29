import uuid as uuid_lib
from datetime import datetime

from django.core.management.base import BaseCommand, CommandError
from django.utils import timezone

from analysis.models import HybridReferenceTask


DEMO_UUID = "3198c38c-6a01-412c-b12c-0610ee1a06b0"
DEMO_TASK_NAME = "demo_task"
DEMO_TCGA_TYPE = "TCGA_ACC"
DEMO_LNCRNA_TYPE = "log2tpm"
DEMO_USE_PADJ = False
DEMO_MAP_INFO = "ImmiRImmiR_ACC"
DEMO_DEG_METHOD = "limma"

DEMO_CREATE_TIME = "2026-06-27 10:27:11"
DEMO_FINISH_TIME = "2026-06-27 10:41:52"

DEMO_MRNA_FILE = "mrna.csv"
DEMO_META_FILE = "meta.csv"

DEMO_LOGFC_CUTOFF_MRNA = 1e-6
DEMO_PADJ_CUTOFF_MRNA = 0.3


def make_aware_datetime(datetime_string: str):
    naive_dt = datetime.strptime(datetime_string, "%Y-%m-%d %H:%M:%S")

    if timezone.is_aware(naive_dt):
        return naive_dt

    return timezone.make_aware(
        naive_dt,
        timezone.get_current_timezone(),
    )


class Command(BaseCommand):
    help = "Create or update a demo HybridReferenceTask record for development."

    def add_arguments(self, parser):
        parser.add_argument(
            "--uuid",
            default=DEMO_UUID,
            help="UUID for the demo task.",
        )

        parser.add_argument(
            "--task-name",
            default=DEMO_TASK_NAME,
            help="Task name for the demo task.",
        )

        parser.add_argument(
            "--user",
            default="",
            help="Optional user field.",
        )

        parser.add_argument(
            "--status",
            type=int,
            choices=[
                HybridReferenceTask.Status.Pending,
                HybridReferenceTask.Status.Running,
                HybridReferenceTask.Status.Success,
                HybridReferenceTask.Status.Failed,
            ],
            default=HybridReferenceTask.Status.Success,
            help="Task status. Default is Success.",
        )

        parser.add_argument(
            "--tcga-type",
            default=DEMO_TCGA_TYPE,
            help="TCGA reference type, for example TCGA_ACC.",
        )

        parser.add_argument(
            "--lncrna-type",
            choices=[
                HybridReferenceTask.LncRNAType.log2count,
                HybridReferenceTask.LncRNAType.log2fpkm,
                HybridReferenceTask.LncRNAType.log2fpkmuq,
                HybridReferenceTask.LncRNAType.log2tpm,
            ],
            default=DEMO_LNCRNA_TYPE,
            help="TCGA lncRNA expression type.",
        )

        parser.add_argument(
            "--use-padj",
            dest="use_padj",
            action="store_true",
            help="Use adjusted P-value threshold.",
        )

        parser.add_argument(
            "--no-use-padj",
            dest="use_padj",
            action="store_false",
            help="Use raw P-value threshold instead of adjusted P-value.",
        )

        parser.set_defaults(use_padj=DEMO_USE_PADJ)

        parser.add_argument(
            "--map-info",
            default=DEMO_MAP_INFO,
            help="Immune annotation map_info value.",
        )

        parser.add_argument(
            "--deg-method",
            choices=[
                HybridReferenceTask.DEGMethod.limma,
                HybridReferenceTask.DEGMethod.deseq2,
            ],
            default=DEMO_DEG_METHOD,
            help="DEG method. Default is limma.",
        )

        parser.add_argument(
            "--mrna-file",
            default=DEMO_MRNA_FILE,
            help="Stored mRNA expression filename.",
        )

        parser.add_argument(
            "--meta-file",
            default=DEMO_META_FILE,
            help="Stored metadata filename.",
        )

        parser.add_argument(
            "--logfc-cutoff-mrna",
            type=float,
            default=DEMO_LOGFC_CUTOFF_MRNA,
            help="mRNA logFC cutoff.",
        )

        parser.add_argument(
            "--padj-cutoff-mrna",
            type=float,
            default=DEMO_PADJ_CUTOFF_MRNA,
            help="mRNA adjusted P-value cutoff.",
        )

        parser.add_argument(
            "--create-time",
            default=DEMO_CREATE_TIME,
            help="Create time, format: YYYY-MM-DD HH:MM:SS.",
        )

        parser.add_argument(
            "--finish-time",
            default=DEMO_FINISH_TIME,
            help="Finish time, format: YYYY-MM-DD HH:MM:SS.",
        )

    def handle(self, *args, **options):
        try:
            task_uuid = uuid_lib.UUID(options["uuid"])
        except ValueError as e:
            raise CommandError(f"Invalid UUID: {options['uuid']}") from e

        task_status = options["status"]
        create_time = make_aware_datetime(options["create_time"])

        if task_status in {
            HybridReferenceTask.Status.Success,
            HybridReferenceTask.Status.Failed,
        }:
            finish_time = make_aware_datetime(options["finish_time"])
        else:
            finish_time = None

        task, created = HybridReferenceTask.objects.update_or_create(
            uuid=task_uuid,
            defaults={
                "user": options["user"],
                "task_name": options["task_name"],
                "status": task_status,
                "tcga_type": options["tcga_type"],
                "lncrna_type": options["lncrna_type"],
                "use_padj": options["use_padj"],
                "map_info": options["map_info"],
                "deg_method": options["deg_method"],
                "mrna_file": options["mrna_file"],
                "meta_file": options["meta_file"],
                "logfc_cutoff_mrna": options["logfc_cutoff_mrna"],
                "padj_cutoff_mrna": options["padj_cutoff_mrna"],
                "finish_time": finish_time,
            },
        )

        HybridReferenceTask.objects.filter(uuid=task_uuid).update(
            create_time=create_time,
        )

        task.refresh_from_db()

        action = "created" if created else "updated"

        self.stdout.write(
            self.style.SUCCESS(
                f"Demo HybridReferenceTask {action}: {task.uuid}"
            )
        )

        self.stdout.write(f"task_name: {task.task_name}")
        self.stdout.write(f"user: {task.user}")
        self.stdout.write(f"status: {task.status} ({task.get_status_display()})")
        self.stdout.write(f"create_time: {task.create_time}")
        self.stdout.write(f"finish_time: {task.finish_time}")

        self.stdout.write(f"tcga_type: {task.tcga_type}")
        self.stdout.write(f"lncrna_type: {task.lncrna_type}")
        self.stdout.write(f"use_padj: {task.use_padj}")
        self.stdout.write(f"map_info: {task.map_info}")
        self.stdout.write(f"deg_method: {task.deg_method}")

        self.stdout.write(f"mrna_file: {task.mrna_file}")
        self.stdout.write(f"meta_file: {task.meta_file}")

        self.stdout.write(f"logfc_cutoff_mrna: {task.logfc_cutoff_mrna}")
        self.stdout.write(f"padj_cutoff_mrna: {task.padj_cutoff_mrna}")

        self.stdout.write(f"workspace: {task.get_workspace_dir_absolute_path()}")
        self.stdout.write(f"input_dir: {task.get_input_dir_absolute_path()}")
        self.stdout.write(f"output_dir: {task.get_output_dir_absolute_path()}")
