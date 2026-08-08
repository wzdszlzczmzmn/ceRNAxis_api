import uuid as uuid_lib

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from django.utils import timezone
from django.utils.dateparse import parse_datetime

from analysis.models import SCSTHybridReferenceTask


DEMO_UUID = "65f4f05b-1fd6-4655-8a93-c1011a513a2c"
DEMO_USER = ""
DEMO_TASK_NAME = "demo_task"
DEMO_STATUS = SCSTHybridReferenceTask.Status.Success

DEMO_DATA_TYPE = SCSTHybridReferenceTask.DataType.SC
DEMO_TCGA_TYPE = "TCGA_BRCA"
DEMO_LNCRNA_TYPE = SCSTHybridReferenceTask.LncRNAType.log2tpm

DEMO_EXP_FILE = "expression.h5ad"
DEMO_META_FILE = ""
DEMO_GROUP_COL = "Celltype (malignancy)"

DEMO_MAP_INFO = "ImmiRImmiR_BRCA"
DEMO_USE_PADJ = False

DEMO_LOGFC_CUTOFF_MRNA = 1e-6
DEMO_PADJ_CUTOFF_MRNA = 0.3

DEMO_CREATE_TIME = "2026-08-08T15:25:58.246846+00:00"
DEMO_FINISH_TIME = "2026-08-08T16:20:22+00:00"


def parse_task_datetime(
    datetime_string: str,
    field_name: str,
):
    """
    Parse an ISO 8601 datetime string and return a timezone-aware datetime.

    Supported examples:
        2026-07-23T14:33:31.464Z
        2026-07-23T15:21:38Z
        2026-07-23T15:21:38+00:00
        2026-07-23 23:21:38
    """
    datetime_string = str(
        datetime_string or ""
    ).strip()

    if not datetime_string:
        raise CommandError(
            f"Missing datetime value: {field_name}."
        )

    parsed_datetime = parse_datetime(
        datetime_string
    )

    if parsed_datetime is None:
        raise CommandError(
            f"Invalid {field_name} value: "
            f"{datetime_string}. "
            "Use an ISO 8601 datetime such as "
            "'2026-07-23T14:33:31.464Z'."
        )

    if timezone.is_naive(parsed_datetime):
        parsed_datetime = timezone.make_aware(
            parsed_datetime,
            timezone.get_current_timezone(),
        )

    return parsed_datetime


class Command(BaseCommand):
    help = (
        "Create or update a demo SC/ST Hybrid Reference task "
        "for development."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--uuid",
            default=DEMO_UUID,
            help="UUID for the demo task.",
        )

        parser.add_argument(
            "--user",
            default=DEMO_USER,
            help="Optional user field.",
        )

        parser.add_argument(
            "--task-name",
            default=DEMO_TASK_NAME,
            help="Task name for the demo task.",
        )

        parser.add_argument(
            "--status",
            type=int,
            choices=[
                SCSTHybridReferenceTask.Status.Pending,
                SCSTHybridReferenceTask.Status.Running,
                SCSTHybridReferenceTask.Status.Success,
                SCSTHybridReferenceTask.Status.Failed,
            ],
            default=DEMO_STATUS,
            help=(
                "Task status: "
                "0=Pending, 1=Running, 2=Success, 3=Failed. "
                "Default is Success."
            ),
        )

        parser.add_argument(
            "--data-type",
            choices=[
                SCSTHybridReferenceTask.DataType.SC,
                SCSTHybridReferenceTask.DataType.ST,
            ],
            default=DEMO_DATA_TYPE,
            help=(
                "Input data type: "
                "'sc' for single-cell RNA-seq or "
                "'st' for spatial transcriptomics."
            ),
        )

        parser.add_argument(
            "--tcga-type",
            default=DEMO_TCGA_TYPE,
            help=(
                "TCGA reference type, for example "
                "TCGA_BRCA."
            ),
        )

        parser.add_argument(
            "--lncrna-type",
            choices=[
                SCSTHybridReferenceTask.LncRNAType.log2count,
                SCSTHybridReferenceTask.LncRNAType.log2fpkm,
                SCSTHybridReferenceTask.LncRNAType.log2fpkmuq,
                SCSTHybridReferenceTask.LncRNAType.log2tpm,
            ],
            default=DEMO_LNCRNA_TYPE,
            help="TCGA lncRNA expression type.",
        )

        parser.add_argument(
            "--exp-file",
            default=DEMO_EXP_FILE,
            help="Stored SC/ST AnnData H5AD filename.",
        )

        parser.add_argument(
            "--meta-file",
            default=DEMO_META_FILE,
            help="Stored metadata filename. Empty for H5AD-only SC/ST input.",
        )

        parser.add_argument(
            "--group-col",
            default=DEMO_GROUP_COL,
            help=(
                "adata.obs column used to group cells "
                "or spatial observations."
            ),
        )

        parser.add_argument(
            "--map-info",
            default=DEMO_MAP_INFO,
            help="Immune annotation map_info value.",
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
            help=(
                "Use raw P-value threshold instead of "
                "adjusted P-value."
            ),
        )

        parser.set_defaults(
            use_padj=DEMO_USE_PADJ
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
            help=(
                "Task creation time in ISO 8601 format. "
                "Example: 2026-07-23T14:33:31.464Z."
            ),
        )

        parser.add_argument(
            "--finish-time",
            default=DEMO_FINISH_TIME,
            help=(
                "Task finish time in ISO 8601 format. "
                "Example: 2026-07-23T15:21:38Z."
            ),
        )

    @transaction.atomic
    def handle(self, *args, **options):
        try:
            task_uuid = uuid_lib.UUID(
                str(options["uuid"])
            )
        except (
            ValueError,
            AttributeError,
            TypeError,
        ) as exc:
            raise CommandError(
                f"Invalid UUID: {options['uuid']}"
            ) from exc

        task_status = options["status"]

        create_time = parse_task_datetime(
            options["create_time"],
            field_name="create_time",
        )

        if task_status in {
            SCSTHybridReferenceTask.Status.Success,
            SCSTHybridReferenceTask.Status.Failed,
        }:
            finish_time = parse_task_datetime(
                options["finish_time"],
                field_name="finish_time",
            )
        else:
            finish_time = None

        task, created = (
            SCSTHybridReferenceTask.objects.update_or_create(
                uuid=task_uuid,
                defaults={
                    "user": options["user"],
                    "task_name": options["task_name"],
                    "status": task_status,
                    "data_type": options["data_type"],
                    "tcga_type": options["tcga_type"],
                    "lncrna_type": options["lncrna_type"],
                    "exp_file": options["exp_file"],
                    "meta_file": options["meta_file"],
                    "group_col": options["group_col"],
                    "map_info": options["map_info"],
                    "use_padj": options["use_padj"],
                    "logfc_cutoff_mrna": options[
                        "logfc_cutoff_mrna"
                    ],
                    "padj_cutoff_mrna": options[
                        "padj_cutoff_mrna"
                    ],
                    "finish_time": finish_time,
                },
            )
        )

        # create_time 使用 auto_now_add=True。
        # save() 时 Django 会自动设置该字段，因此需要通过
        # QuerySet.update() 覆盖为原始 demo 任务时间。
        SCSTHybridReferenceTask.objects.filter(
            pk=task.pk,
        ).update(
            create_time=create_time,
        )

        task.refresh_from_db()

        action = (
            "created"
            if created
            else "updated"
        )

        self.stdout.write(
            self.style.SUCCESS(
                "Demo SCSTHybridReferenceTask "
                f"{action}: {task.uuid}"
            )
        )

        self.stdout.write(
            f"id: {task.id}"
        )
        self.stdout.write(
            f"uuid: {task.uuid}"
        )
        self.stdout.write(
            f"user: {task.user}"
        )
        self.stdout.write(
            f"task_name: {task.task_name}"
        )
        self.stdout.write(
            f"status: {task.status} "
            f"({task.get_status_display()})"
        )

        self.stdout.write(
            f"data_type: {task.data_type} "
            f"({task.get_data_type_display()})"
        )
        self.stdout.write(
            f"tcga_type: {task.tcga_type}"
        )
        self.stdout.write(
            f"lncrna_type: {task.lncrna_type}"
        )

        self.stdout.write(
            f"exp_file: {task.exp_file}"
        )
        self.stdout.write(
            f"meta_file: {task.meta_file}"
        )
        self.stdout.write(
            f"group_col: {task.group_col}"
        )

        self.stdout.write(
            f"map_info: {task.map_info}"
        )
        self.stdout.write(
            f"use_padj: {task.use_padj}"
        )

        self.stdout.write(
            "logfc_cutoff_mrna: "
            f"{task.logfc_cutoff_mrna}"
        )
        self.stdout.write(
            "padj_cutoff_mrna: "
            f"{task.padj_cutoff_mrna}"
        )

        self.stdout.write(
            f"create_time: {task.create_time}"
        )
        self.stdout.write(
            f"finish_time: {task.finish_time}"
        )

        self.stdout.write(
            "workspace: "
            f"{task.get_workspace_dir_absolute_path()}"
        )
        self.stdout.write(
            "input_dir: "
            f"{task.get_input_dir_absolute_path()}"
        )
        self.stdout.write(
            "output_dir: "
            f"{task.get_output_dir_absolute_path()}"
        )
