import uuid as uuid_lib
from datetime import datetime

from django.core.management.base import BaseCommand, CommandError
from django.utils import timezone

from analysis.models import PairedCohortTask


DEMO_UUID = "9de2c6b5-034c-4c4c-811a-976fbac53215"
DEMO_TASK_NAME = "demo_task"
DEMO_MAP_INFO = "ImmiRImmiR_ACC"
DEMO_DEG_METHOD = "limma"
DEMO_CREATE_TIME = "2026-06-16 09:24:09"
DEMO_FINISH_TIME = "2026-06-16 09:24:48"

DEMO_MRNA_FILE = "mrna.csv"
DEMO_MIRNA_FILE = "mirna.csv"
DEMO_LNCRNA_FILE = "lncrna.csv"
DEMO_META_FILE = "meta.csv"

DEMO_LOGFC_CUTOFF_MRNA = 1.0
DEMO_PADJ_CUTOFF_MRNA = 0.1

DEMO_LOGFC_CUTOFF_MIRNA = 0.5
DEMO_PADJ_CUTOFF_MIRNA = 0.3

DEMO_LOGFC_CUTOFF_LNCRNA = 0.5
DEMO_PADJ_CUTOFF_LNCRNA = 0.3


def make_aware_datetime(datetime_string: str):
    naive_dt = datetime.strptime(datetime_string, "%Y-%m-%d %H:%M:%S")

    if timezone.is_aware(naive_dt):
        return naive_dt

    return timezone.make_aware(
        naive_dt,
        timezone.get_current_timezone(),
    )


class Command(BaseCommand):
    help = "Create or update a demo PairedCohortTask record for development."

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
            "--map-info",
            default=DEMO_MAP_INFO,
            help="Immune annotation map_info value.",
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
                PairedCohortTask.Status.Success,
                PairedCohortTask.Status.Pending,
                PairedCohortTask.Status.Failed,
                PairedCohortTask.Status.Running,
            ],
            default=PairedCohortTask.Status.Success,
            help="Task status. Default is Success.",
        )

        parser.add_argument(
            "--deg-method",
            choices=["limma", "deseq2"],
            default=DEMO_DEG_METHOD,
            help="DEG method. Default is limma.",
        )

        parser.add_argument(
            "--mrna-file",
            default=DEMO_MRNA_FILE,
            help="Stored mRNA expression filename.",
        )

        parser.add_argument(
            "--mirna-file",
            default=DEMO_MIRNA_FILE,
            help="Stored miRNA expression filename.",
        )

        parser.add_argument(
            "--lncrna-file",
            default=DEMO_LNCRNA_FILE,
            help="Stored lncRNA expression filename.",
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
            "--logfc-cutoff-mirna",
            type=float,
            default=DEMO_LOGFC_CUTOFF_MIRNA,
            help="miRNA logFC cutoff.",
        )

        parser.add_argument(
            "--padj-cutoff-mirna",
            type=float,
            default=DEMO_PADJ_CUTOFF_MIRNA,
            help="miRNA adjusted P-value cutoff.",
        )

        parser.add_argument(
            "--logfc-cutoff-lncrna",
            type=float,
            default=DEMO_LOGFC_CUTOFF_LNCRNA,
            help="lncRNA logFC cutoff.",
        )

        parser.add_argument(
            "--padj-cutoff-lncrna",
            type=float,
            default=DEMO_PADJ_CUTOFF_LNCRNA,
            help="lncRNA adjusted P-value cutoff.",
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
            PairedCohortTask.Status.Success,
            PairedCohortTask.Status.Failed,
        }:
            finish_time = make_aware_datetime(options["finish_time"])
        else:
            finish_time = None

        task, created = PairedCohortTask.objects.update_or_create(
            uuid=task_uuid,
            defaults={
                "user": options["user"],
                "task_name": options["task_name"],
                "status": task_status,
                "finish_time": finish_time,
                "map_info": options["map_info"],
                "deg_method": options["deg_method"],
                "mrna_file": options["mrna_file"],
                "mirna_file": options["mirna_file"],
                "lncrna_file": options["lncrna_file"],
                "meta_file": options["meta_file"],
                "logfc_cutoff_mrna": options["logfc_cutoff_mrna"],
                "padj_cutoff_mrna": options["padj_cutoff_mrna"],
                "logfc_cutoff_mirna": options["logfc_cutoff_mirna"],
                "padj_cutoff_mirna": options["padj_cutoff_mirna"],
                "logfc_cutoff_lncrna": options["logfc_cutoff_lncrna"],
                "padj_cutoff_lncrna": options["padj_cutoff_lncrna"],
            },
        )

        PairedCohortTask.objects.filter(uuid=task_uuid).update(
            create_time=create_time,
        )

        task.refresh_from_db()

        action = "created" if created else "updated"

        self.stdout.write(
            self.style.SUCCESS(
                f"Demo PairedCohortTask {action}: {task.uuid}"
            )
        )
        self.stdout.write(f"task_name: {task.task_name}")
        self.stdout.write(f"user: {task.user}")
        self.stdout.write(f"status: {task.status} ({task.get_status_display()})")
        self.stdout.write(f"create_time: {task.create_time}")
        self.stdout.write(f"finish_time: {task.finish_time}")
        self.stdout.write(f"map_info: {task.map_info}")
        self.stdout.write(f"deg_method: {task.deg_method}")
        self.stdout.write(f"mrna_file: {task.mrna_file}")
        self.stdout.write(f"mirna_file: {task.mirna_file}")
        self.stdout.write(f"lncrna_file: {task.lncrna_file}")
        self.stdout.write(f"meta_file: {task.meta_file}")
        self.stdout.write(f"logfc_cutoff_mrna: {task.logfc_cutoff_mrna}")
        self.stdout.write(f"padj_cutoff_mrna: {task.padj_cutoff_mrna}")
        self.stdout.write(f"logfc_cutoff_mirna: {task.logfc_cutoff_mirna}")
        self.stdout.write(f"padj_cutoff_mirna: {task.padj_cutoff_mirna}")
        self.stdout.write(f"logfc_cutoff_lncrna: {task.logfc_cutoff_lncrna}")
        self.stdout.write(f"padj_cutoff_lncrna: {task.padj_cutoff_lncrna}")
        self.stdout.write(f"workspace: {task.get_workspace_dir_absolute_path()}")
        self.stdout.write(f"input_dir: {task.get_input_dir_absolute_path()}")
        self.stdout.write(f"output_dir: {task.get_output_dir_absolute_path()}")
