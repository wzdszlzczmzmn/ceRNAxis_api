import uuid as uuid_lib

from django.core.management.base import BaseCommand, CommandError
from django.utils import timezone
from django.utils.dateparse import parse_datetime

from analysis.models import PairedCohortTask


DEMO_UUID = "37051bf6-5b88-4ec8-a958-b36a28070109"
DEMO_TASK_NAME = "demo_task"
DEMO_MAP_INFO = "ImmiRImmiR_LUAD"
DEMO_DEG_METHOD = "limma"
DEMO_CANCER_TYPE = "LUAD"
DEMO_USE_PADJ = False

DEMO_CREATE_TIME = "2026-07-22T03:20:39.487Z"
DEMO_FINISH_TIME = "2026-07-22T03:40:03Z"

DEMO_MRNA_FILE = "mrna.csv"
DEMO_MIRNA_FILE = "mirna.csv"
DEMO_LNCRNA_FILE = "lncrna.csv"
DEMO_CIRCRNA_FILE = "circrna.csv"
DEMO_META_FILE = "meta.csv"

DEMO_LOGFC_CUTOFF_MRNA = 1e-6
DEMO_PADJ_CUTOFF_MRNA = 0.1

DEMO_LOGFC_CUTOFF_MIRNA = 1e-6
DEMO_PADJ_CUTOFF_MIRNA = 0.3

DEMO_LOGFC_CUTOFF_LNCRNA = 1e-6
DEMO_PADJ_CUTOFF_LNCRNA = 0.3

DEMO_LOGFC_CUTOFF_CIRCRNA = 1e-6
DEMO_PADJ_CUTOFF_CIRCRNA = 0.3


def parse_task_datetime(datetime_string: str):
    """
    Parse an ISO 8601 datetime string and return a timezone-aware datetime.

    Supported examples:
        2026-07-22T03:20:39.487Z
        2026-07-22T03:40:03Z
        2026-07-22T03:40:03+00:00
        2026-07-22 12:40:03
    """
    parsed_datetime = parse_datetime(datetime_string)

    if parsed_datetime is None:
        raise CommandError(
            f"Invalid datetime value: {datetime_string}. "
            "Use an ISO 8601 datetime such as "
            "'2026-07-22T03:20:39.487Z'."
        )

    if timezone.is_naive(parsed_datetime):
        parsed_datetime = timezone.make_aware(
            parsed_datetime,
            timezone.get_current_timezone(),
        )

    return parsed_datetime


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
                PairedCohortTask.Status.Pending,
                PairedCohortTask.Status.Running,
                PairedCohortTask.Status.Success,
                PairedCohortTask.Status.Failed,
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
            "--cancer-type",
            default=DEMO_CANCER_TYPE,
            help="Cancer type for the demo task.",
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
            "--circrna-file",
            default=DEMO_CIRCRNA_FILE,
            help="Stored circRNA expression filename.",
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
            "--logfc-cutoff-circrna",
            type=float,
            default=DEMO_LOGFC_CUTOFF_CIRCRNA,
            help="circRNA logFC cutoff.",
        )

        parser.add_argument(
            "--padj-cutoff-circrna",
            type=float,
            default=DEMO_PADJ_CUTOFF_CIRCRNA,
            help="circRNA adjusted P-value cutoff.",
        )

        parser.add_argument(
            "--create-time",
            default=DEMO_CREATE_TIME,
            help=(
                "Task creation time in ISO 8601 format. "
                "Example: 2026-07-22T03:20:39.487Z."
            ),
        )

        parser.add_argument(
            "--finish-time",
            default=DEMO_FINISH_TIME,
            help=(
                "Task finish time in ISO 8601 format. "
                "Example: 2026-07-22T03:40:03Z."
            ),
        )

    def handle(self, *args, **options):
        try:
            task_uuid = uuid_lib.UUID(options["uuid"])
        except (ValueError, AttributeError, TypeError) as exc:
            raise CommandError(
                f"Invalid UUID: {options['uuid']}"
            ) from exc

        task_status = options["status"]
        create_time = parse_task_datetime(options["create_time"])

        if task_status in {
            PairedCohortTask.Status.Success,
            PairedCohortTask.Status.Failed,
        }:
            finish_time = parse_task_datetime(options["finish_time"])
        else:
            finish_time = None

        task, created = PairedCohortTask.objects.update_or_create(
            uuid=task_uuid,
            defaults={
                "user": options["user"],
                "task_name": options["task_name"],
                "status": task_status,
                "map_info": options["map_info"],
                "deg_method": options["deg_method"],
                "cancer_type": options["cancer_type"],
                "use_padj": options["use_padj"],
                "mrna_file": options["mrna_file"],
                "mirna_file": options["mirna_file"],
                "lncrna_file": options["lncrna_file"],
                "circrna_file": options["circrna_file"],
                "meta_file": options["meta_file"],
                "logfc_cutoff_mrna": options["logfc_cutoff_mrna"],
                "padj_cutoff_mrna": options["padj_cutoff_mrna"],
                "logfc_cutoff_mirna": options["logfc_cutoff_mirna"],
                "padj_cutoff_mirna": options["padj_cutoff_mirna"],
                "logfc_cutoff_lncrna": options["logfc_cutoff_lncrna"],
                "padj_cutoff_lncrna": options["padj_cutoff_lncrna"],
                "logfc_cutoff_circrna": options[
                    "logfc_cutoff_circrna"
                ],
                "padj_cutoff_circrna": options[
                    "padj_cutoff_circrna"
                ],
                "finish_time": finish_time,
            },
        )

        # create_time 使用 auto_now_add=True，不能通过 save() 正常覆盖。
        # 因此在记录创建或更新后，通过 QuerySet.update() 设置原始时间。
        PairedCohortTask.objects.filter(
            uuid=task_uuid,
        ).update(
            create_time=create_time,
        )

        task.refresh_from_db()

        action = "created" if created else "updated"

        self.stdout.write(
            self.style.SUCCESS(
                f"Demo PairedCohortTask {action}: {task.uuid}"
            )
        )

        self.stdout.write(f"id: {task.id}")
        self.stdout.write(f"uuid: {task.uuid}")
        self.stdout.write(f"user: {task.user}")
        self.stdout.write(f"task_name: {task.task_name}")
        self.stdout.write(
            f"status: {task.status} ({task.get_status_display()})"
        )
        self.stdout.write(f"map_info: {task.map_info}")
        self.stdout.write(f"deg_method: {task.deg_method}")
        self.stdout.write(f"cancer_type: {task.cancer_type}")
        self.stdout.write(f"use_padj: {task.use_padj}")

        self.stdout.write(f"mrna_file: {task.mrna_file}")
        self.stdout.write(f"mirna_file: {task.mirna_file}")
        self.stdout.write(f"lncrna_file: {task.lncrna_file}")
        self.stdout.write(f"circrna_file: {task.circrna_file}")
        self.stdout.write(f"meta_file: {task.meta_file}")

        self.stdout.write(
            f"logfc_cutoff_mrna: {task.logfc_cutoff_mrna}"
        )
        self.stdout.write(
            f"padj_cutoff_mrna: {task.padj_cutoff_mrna}"
        )

        self.stdout.write(
            f"logfc_cutoff_mirna: {task.logfc_cutoff_mirna}"
        )
        self.stdout.write(
            f"padj_cutoff_mirna: {task.padj_cutoff_mirna}"
        )

        self.stdout.write(
            f"logfc_cutoff_lncrna: {task.logfc_cutoff_lncrna}"
        )
        self.stdout.write(
            f"padj_cutoff_lncrna: {task.padj_cutoff_lncrna}"
        )

        self.stdout.write(
            f"logfc_cutoff_circrna: {task.logfc_cutoff_circrna}"
        )
        self.stdout.write(
            f"padj_cutoff_circrna: {task.padj_cutoff_circrna}"
        )

        self.stdout.write(f"create_time: {task.create_time}")
        self.stdout.write(f"finish_time: {task.finish_time}")

        self.stdout.write(
            f"workspace: {task.get_workspace_dir_absolute_path()}"
        )
        self.stdout.write(
            f"input_dir: {task.get_input_dir_absolute_path()}"
        )
        self.stdout.write(
            f"output_dir: {task.get_output_dir_absolute_path()}"
        )
