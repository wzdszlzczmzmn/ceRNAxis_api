from database.models import (
    AxisDatasetSource,
    AxisGroupType,
    AxisModule,
)

from .import_contracts import AxisContextSpec
from .import_normalization import (
    AxisImportValidationError,
)


def validate_context_spec(
    spec: AxisContextSpec,
) -> None:
    if spec.module == AxisModule.MODULE2:
        if (
            spec.dataset_source
            != AxisDatasetSource.TCGA
        ):
            raise AxisImportValidationError(
                "Module 2 context must use "
                "dataset_source='TCGA'."
            )

        if not spec.dataset_name.endswith("_mRNA"):
            raise AxisImportValidationError(
                "Module 2 dataset_name must end "
                "with '_mRNA'."
            )

        if (
            spec.group_type != AxisGroupType.NONE
            or spec.group_by
        ):
            raise AxisImportValidationError(
                "Module 2 context must use "
                "group_type='none' and empty group_by."
            )

        return

    if spec.module == AxisModule.MODULE3:
        if (
            spec.dataset_source
            != AxisDatasetSource.TIMEDB
        ):
            raise AxisImportValidationError(
                "Module 3 context must use "
                "dataset_source='TIMEDB'."
            )

        if spec.group_type == AxisGroupType.NONE:
            raise AxisImportValidationError(
                "Module 3 context must use "
                "other, grade, or stage."
            )

        if not spec.group_by:
            raise AxisImportValidationError(
                "Module 3 context must provide group_by."
            )

        return

    raise AxisImportValidationError(
        f"Unsupported Axis module: {spec.module!r}."
    )
