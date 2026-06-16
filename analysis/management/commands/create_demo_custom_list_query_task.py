import uuid as uuid_lib
from datetime import datetime

from django.core.management.base import BaseCommand, CommandError
from django.utils import timezone

from analysis.models import CustomListQueryTask


DEMO_UUID = "9d8abd76-23b6-40f9-8e7f-cbeb8f11cc83"
DEMO_TASK_NAME = "demo_task"
DEMO_MAP_INFO = "ImmiRImmiR_ACC"
DEMO_FINISH_TIME = "2026-06-12 16:43:07"


DEMO_RNAS = {
    "miRNA": [
        "hsa-mir-567",
        "hsa-miR-6742-5p",
        "hsa-miR-127-5p",
        "hsa-miR-9-3p",
        "hsa-mir-1305",
        "hsa-miR-4487",
        "hsa-mir-629",
        "hsa-mir-181a-2",
        "hsa-miR-548o-3p",
        "hsa-miR-642a-5p",
        "hsa-miR-12115",
        "hsa-miR-5011-5p",
        "hsa-miR-548az-5p",
        "mmu-miR-140-5p",
        "hsa-miR-6715b-5p",
        "hsa-miR-660-3p",
        "hsa-miR-3617-3p",
        "hsa-miR-1266-3p",
        "hsa-miR-639",
        "hsa-miR-10525-3p",
        "hsa-miR-598-3p",
        "hsa-miR-4483",
        "hsa-miR-4694-3p",
        "hsa-miR-548p",
        "hsa-miR-365b-3p",
        "hsa-mir-513c",
        "hsa-mir-4423",
        "hsa-miR-1299",
        "hsa-mir-320a",
        "hsa-mir-421",
    ],
    "mRNA": [
        "LRRC24",
        "CTSC",
        "UBE2K",
        "NDUFB1",
        "PCTP",
        "TLCD5",
        "FHL3",
        "GALNT16",
        "IGF2R",
        "VILL",
        "CDT1",
        "NRARP",
        "APMAP",
        "ZNF532",
        "SYNM",
        "SMG7",
        "KDM8",
        "HLA-DMA",
        "TNP2",
        "ST3GAL4",
        "ZNF544",
        "ACSL6",
        "RRP1B",
        "STX4",
        "HPX",
        "C2CD2L",
        "OSCP1",
        "PAFAH1B3",
        "GALNS",
        "SPDYE14",
    ],
    "lncRNA": [
        "NONHSAG043011",
        "AC092279.1",
        "AC139099.1",
        "RP11-108M9.3",
        "AL645939.2",
        "AC022075.2",
        "NONHSAG054710",
        "LINC01630",
        "AC093515.1",
        "NONHSAG043672",
        "AC113383.1",
        "RP11-274H2.5",
        "RP1-261D10.2",
        "AC021242.3",
        "AL137782.1",
        "LOC440461",
        "GS1-124K5.4",
        "RP11-539L10.3",
        "NONHSAG011425",
        "NONHSAG026847",
        "AP001055.6",
        "RP11-56D16.8",
        "LINC01296",
        "MIR22HG",
        "Z95704.3",
        "NONHSAG036097",
        "AC090607.5",
        "AC016907.3",
        "SNAP47-IT1",
        "AC005220.1",
    ],
    "circRNA": [],
}


def make_aware_datetime(datetime_string: str):
    naive_dt = datetime.strptime(datetime_string, "%Y-%m-%d %H:%M:%S")

    if timezone.is_aware(naive_dt):
        return naive_dt

    return timezone.make_aware(
        naive_dt,
        timezone.get_current_timezone(),
    )


class Command(BaseCommand):
    help = "Create or update a demo CustomListQueryTask record for development."

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
            choices=[
                CustomListQueryTask.Status.Success,
                CustomListQueryTask.Status.Pending,
                CustomListQueryTask.Status.Failed,
                CustomListQueryTask.Status.Running,
            ],
            default=CustomListQueryTask.Status.Success,
            help="Task status. Default is Success.",
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

        if task_status in {
            CustomListQueryTask.Status.Success,
            CustomListQueryTask.Status.Failed,
        }:
            finish_time = make_aware_datetime(options["finish_time"])
        else:
            finish_time = None

        task, created = CustomListQueryTask.objects.update_or_create(
            uuid=task_uuid,
            defaults={
                "user": options["user"],
                "task_name": options["task_name"],
                "status": task_status,
                "finish_time": finish_time,
                "map_info": options["map_info"],
                "rnas": DEMO_RNAS,
            },
        )

        action = "created" if created else "updated"

        self.stdout.write(
            self.style.SUCCESS(
                f"Demo CustomListQueryTask {action}: {task.uuid}"
            )
        )
        self.stdout.write(f"task_name: {task.task_name}")
        self.stdout.write(f"status: {task.status} ({task.get_status_display()})")
        self.stdout.write(f"create_time: {task.create_time}")
        self.stdout.write(f"finish_time: {task.finish_time}")
        self.stdout.write(f"map_info: {task.map_info}")
        self.stdout.write(f"miRNA_count: {task.miRNA_count}")
        self.stdout.write(f"mRNA_count: {task.mRNA_count}")
        self.stdout.write(f"lncRNA_count: {task.lncRNA_count}")
        self.stdout.write(f"circRNA_count: {task.circRNA_count}")
        self.stdout.write(f"total_rna_count: {task.total_rna_count}")
