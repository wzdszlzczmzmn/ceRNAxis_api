from database.models import (
    AxisDatasetSource,
    AxisGroupType,
    AxisModule, AxisResultKind,
)

from .import_contracts import AxisContextSpec
from .import_normalization import (
    AxisImportValidationError,
)


def validate_context_spec(
    spec: AxisContextSpec,
) -> None:
    # -------------------------
    # TCGA
    # -------------------------

    if spec.dataset_source == AxisDatasetSource.TCGA:
        if spec.module != AxisModule.MODULE2:
            raise AxisImportValidationError(
                "TCGA context must use module='module2'."
            )

        if not spec.dataset_name.endswith("_mRNA"):
            raise AxisImportValidationError(
                "TCGA dataset_name must end with '_mRNA'."
            )

        if spec.group_type != AxisGroupType.NONE:
            raise AxisImportValidationError(
                "TCGA context must use group_type='none'."
            )

        if spec.group_by:
            raise AxisImportValidationError(
                "TCGA context must use an empty group_by."
            )

        if spec.group_value:
            raise AxisImportValidationError(
                "TCGA context must use an empty group_value."
            )

        return

    # -------------------------
    # TIMEDB
    # -------------------------

    if spec.dataset_source == AxisDatasetSource.TIMEDB:
        if spec.module != AxisModule.MODULE3:
            raise AxisImportValidationError(
                "TIMEDB context must use module='module3'."
            )

        if spec.group_type == AxisGroupType.NONE:
            raise AxisImportValidationError(
                "TIMEDB context must use "
                "group_type='other', 'grade', or 'stage'."
            )

        if not spec.group_by:
            raise AxisImportValidationError(
                "TIMEDB context must provide group_by."
            )

        if spec.group_value:
            raise AxisImportValidationError(
                "TIMEDB context must use an empty group_value."
            )

        return

    # -------------------------
    # SC / ST
    # -------------------------

    if spec.dataset_source in {
        AxisDatasetSource.SC,
        AxisDatasetSource.ST,
    }:
        if spec.module != AxisModule.MODULE3:
            raise AxisImportValidationError(
                f"{spec.dataset_source} context must use "
                "module='module3'."
            )

        if spec.group_type != AxisGroupType.OTHER:
            raise AxisImportValidationError(
                f"{spec.dataset_source} context must use "
                "group_type='other'."
            )

        if not spec.group_by:
            raise AxisImportValidationError(
                f"{spec.dataset_source} context must provide "
                "group_by."
            )

        if not spec.group_value:
            raise AxisImportValidationError(
                f"{spec.dataset_source} context must provide "
                "group_value."
            )

        return

    raise AxisImportValidationError(
        "Unsupported Axis dataset_source: "
        f"{spec.dataset_source!r}."
    )


def validate_context_dataset_metadata(
    *,
    spec: AxisContextSpec,
    dataset_metadata,
) -> None:
    if (
        spec.dataset_source == AxisDatasetSource.SC
        and dataset_metadata.obs_type != "cell"
    ):
        raise AxisImportValidationError(
            "SC context must reference a DatasetMetadata "
            "with obs_type='cell'."
        )

    if (
        spec.dataset_source == AxisDatasetSource.ST
        and dataset_metadata.obs_type != "spot"
    ):
        raise AxisImportValidationError(
            "ST context must reference a DatasetMetadata "
            "with obs_type='spot'."
        )


def validate_context_result_kind(
    *,
    spec: AxisContextSpec,
    result_kind: str,
) -> None:
    if (
        spec.dataset_source
        in {
            AxisDatasetSource.SC,
            AxisDatasetSource.ST,
        }
        and result_kind != AxisResultKind.AXIS_FINAL
    ):
        raise AxisImportValidationError(
            "SC/ST contexts only support "
            "result_kind='axis_final'."
        )
