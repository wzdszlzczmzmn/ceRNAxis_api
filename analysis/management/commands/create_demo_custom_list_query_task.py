import uuid as uuid_lib
from datetime import datetime

from django.core.management.base import BaseCommand, CommandError
from django.utils import timezone

from analysis.models import CustomListQueryTask


DEMO_UUID = "9720aebe-3be4-4262-bd57-cb02d09ce42f"
DEMO_TASK_NAME = "demo_task"
DEMO_STATUS = CustomListQueryTask.Status.Success

DEMO_USER = ""
DEMO_MAP_INFO = ""
DEMO_CANCER_TYPE = "ACC"
DEMO_HAS_MRNA_DIRECTION = True

DEMO_CREATE_TIME = "2026-07-13T07:06:01.475Z"
DEMO_FINISH_TIME = "2026-07-13T08:27:27Z"


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
    "mRNA": [],
    "mRNA_up": [
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
        "CD8A",
        "HLA-DRA",
        "HLA-DPA1",
        "HLA-DQB1",
        "CD74",
        "IFNG",
        "HLA-DRB5",
        "HLA-DPB1",
        "HLA-DQA1",
        "HLA-DRB1",
        "KLRC1",
        "HLA-DQA2",
        "CTSS",
        "KLRD1",
        "HLA-DOA",
        "CIITA",
        "CD8B",
        "KIR2DL4",
        "KLRC3",
        "CD1B",
        "HLA-DMB",
        "CTSB",
        "CD209",
        "CD4",
        "HLA-B",
        "MICB",
        "HLA-C",
        "B2M",
    ],
    "mRNA_down": [
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
        "SH2D1A",
        "IFNG",
        "GZMB",
        "LCP2",
        "PTPN6",
        "PIK3CD",
        "ITGB2",
        "LCK",
        "FASLG",
        "KLRD1",
        "PRF1",
        "HLA-E",
        "CD247",
        "CD48",
        "CD244",
        "KIR2DL4",
        "KLRC3",
        "IFNB1",
        "KLRC1",
        "MICB",
        "SH2D1B",
        "FCER1G",
        "HCST",
        "HLA-B",
        "HLA-C",
        "HLA-A",
        "TYROBP",
        "KLRC2",
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


def parse_datetime(datetime_string: str):
    """
    解析以下格式：

    2026-07-13T07:06:01.475Z
    2026-07-13T07:06:01Z
    2026-07-13T07:06:01+00:00
    2026-07-13 07:06:01
    """
    normalized_value = datetime_string.strip()

    if normalized_value.endswith("Z"):
        normalized_value = f"{normalized_value[:-1]}+00:00"

    try:
        parsed_datetime = datetime.fromisoformat(normalized_value)
    except ValueError as exc:
        raise CommandError(
            f"Invalid datetime value: {datetime_string}"
        ) from exc

    if timezone.is_naive(parsed_datetime):
        parsed_datetime = timezone.make_aware(
            parsed_datetime,
            timezone.get_current_timezone(),
        )

    return parsed_datetime


class Command(BaseCommand):
    help = "Create or update a demo CustomListQueryTask record."

    def add_arguments(self, parser):
        parser.add_argument(
            "--uuid",
            default=DEMO_UUID,
            help="UUID for the task.",
        )

        parser.add_argument(
            "--task-name",
            default=DEMO_TASK_NAME,
            help="Task name.",
        )

        parser.add_argument(
            "--user",
            default=DEMO_USER,
            help="Optional user value.",
        )

        parser.add_argument(
            "--status",
            choices=[
                CustomListQueryTask.Status.Success,
                CustomListQueryTask.Status.Pending,
                CustomListQueryTask.Status.Failed,
                CustomListQueryTask.Status.Running,
            ],
            default=DEMO_STATUS,
            help="Task status.",
        )

        parser.add_argument(
            "--map-info",
            default=DEMO_MAP_INFO,
            help="Deprecated map_info value.",
        )

        parser.add_argument(
            "--cancer-type",
            default=DEMO_CANCER_TYPE,
            help="Cancer type, for example ACC.",
        )

        parser.add_argument(
            "--has-mrna-direction",
            dest="has_mrna_direction",
            action="store_true",
            help="Use directional mRNA input.",
        )

        parser.add_argument(
            "--no-has-mrna-direction",
            dest="has_mrna_direction",
            action="store_false",
            help="Use standard non-directional mRNA input.",
        )

        parser.set_defaults(
            has_mrna_direction=DEMO_HAS_MRNA_DIRECTION,
        )

        parser.add_argument(
            "--create-time",
            default=DEMO_CREATE_TIME,
            help="Task creation time.",
        )

        parser.add_argument(
            "--finish-time",
            default=DEMO_FINISH_TIME,
            help="Task finish time.",
        )

    def handle(self, *args, **options):
        try:
            task_uuid = uuid_lib.UUID(options["uuid"])
        except ValueError as exc:
            raise CommandError(
                f"Invalid UUID: {options['uuid']}"
            ) from exc

        task_status = options["status"]
        create_time = parse_datetime(options["create_time"])

        if task_status in {
            CustomListQueryTask.Status.Success,
            CustomListQueryTask.Status.Failed,
        }:
            finish_time = parse_datetime(options["finish_time"])
        else:
            finish_time = None

        self.validate_rnas(
            rnas=DEMO_RNAS,
            has_mrna_direction=options["has_mrna_direction"],
        )

        task, created = CustomListQueryTask.objects.update_or_create(
            uuid=task_uuid,
            defaults={
                "user": options["user"],
                "task_name": options["task_name"],
                "status": task_status,
                "finish_time": finish_time,
                "map_info": options["map_info"],
                "cancer_type": options["cancer_type"],
                "has_mrna_direction": options["has_mrna_direction"],
                "rnas": DEMO_RNAS,
            },
        )

        # create_time 使用 auto_now_add，不能通过 update_or_create defaults 覆盖，
        # 因此在记录创建或更新后单独执行 update。
        CustomListQueryTask.objects.filter(
            uuid=task_uuid,
        ).update(
            create_time=create_time,
        )

        task.refresh_from_db()

        action = "created" if created else "updated"

        self.stdout.write(
            self.style.SUCCESS(
                f"CustomListQueryTask {action}: {task.uuid}"
            )
        )

        self.stdout.write(f"task_name: {task.task_name}")
        self.stdout.write(f"user: {task.user}")
        self.stdout.write(
            f"status: {task.status} ({task.get_status_display()})"
        )
        self.stdout.write(f"create_time: {task.create_time}")
        self.stdout.write(f"finish_time: {task.finish_time}")
        self.stdout.write(f"map_info: {task.map_info}")
        self.stdout.write(f"cancer_type: {task.cancer_type}")
        self.stdout.write(
            f"has_mrna_direction: {task.has_mrna_direction}"
        )

        self.stdout.write(f"miRNA_count: {task.miRNA_count}")
        self.stdout.write(f"mRNA_count: {task.mRNA_count}")
        self.stdout.write(f"mRNA_up_count: {task.mRNA_up_count}")
        self.stdout.write(f"mRNA_down_count: {task.mRNA_down_count}")
        self.stdout.write(f"lncRNA_count: {task.lncRNA_count}")
        self.stdout.write(f"circRNA_count: {task.circRNA_count}")
        self.stdout.write(f"total_rna_count: {task.total_rna_count}")

        self.stdout.write(
            f"workspace: {task.get_workspace_dir_absolute_path()}"
        )
        self.stdout.write(
            f"input_dir: {task.get_input_dir_absolute_path()}"
        )
        self.stdout.write(
            f"output_dir: {task.get_output_dir_absolute_path()}"
        )

    @staticmethod
    def validate_rnas(rnas: dict, has_mrna_direction: bool):
        required_keys = {
            "miRNA",
            "mRNA",
            "mRNA_up",
            "mRNA_down",
            "lncRNA",
            "circRNA",
        }

        missing_keys = required_keys - set(rnas)

        if missing_keys:
            raise CommandError(
                f"Missing RNA keys: {sorted(missing_keys)}"
            )

        for rna_type in required_keys:
            if not isinstance(rnas[rna_type], list):
                raise CommandError(
                    f"RNA value must be a list: {rna_type}"
                )

        if has_mrna_direction:
            if rnas["mRNA"]:
                raise CommandError(
                    "mRNA must be empty when has_mrna_direction is true."
                )
        else:
            if rnas["mRNA_up"] or rnas["mRNA_down"]:
                raise CommandError(
                    "mRNA_up and mRNA_down must be empty when "
                    "has_mrna_direction is false."
                )
