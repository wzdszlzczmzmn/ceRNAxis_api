from django.db import transaction
from django.db.models import Count, Q

from database.models import (
    AxisContextPresence,
    AxisDatasetContext,
    AxisObservation,
    AxisResultKind,
)


@transaction.atomic
def rebuild_axis_context_presence(
    *,
    context: AxisDatasetContext,
    batch_size: int = 5000,
) -> dict:
    deleted_count, _ = (
        AxisContextPresence.objects
        .filter(context=context)
        .delete()
    )

    if not context.is_active:
        return {
            "context_id": context.id,
            "deleted_count": deleted_count,
            "created_count": 0,
        }

    rows = (
        AxisObservation.objects
        .filter(
            artifact__context=context,
            artifact__is_active=True,
        )
        .values("axis_id")
        .annotate(
            observation_count=Count("id"),

            axis_final_observation_count=Count(
                "id",
                filter=Q(
                    artifact__result_kind=(
                        AxisResultKind.AXIS_FINAL
                    ),
                ),
            ),

            sponge_observation_count=Count(
                "id",
                filter=Q(
                    artifact__result_kind=(
                        AxisResultKind.SPONGE
                    ),
                ),
            ),
        )
        .order_by("axis_id")
    )

    objects = []

    for row in rows.iterator(
        chunk_size=batch_size,
    ):
        axis_final_count = int(
            row["axis_final_observation_count"]
        )
        sponge_count = int(
            row["sponge_observation_count"]
        )

        objects.append(
            AxisContextPresence(
                context=context,
                axis_id=row["axis_id"],
                observation_count=int(
                    row["observation_count"]
                ),
                axis_final_observation_count=(
                    axis_final_count
                ),
                sponge_observation_count=(
                    sponge_count
                ),
                has_axis_final=axis_final_count > 0,
                has_sponge=sponge_count > 0,
            )
        )

    created = AxisContextPresence.objects.bulk_create(
        objects,
        batch_size=batch_size,
    )

    return {
        "context_id": context.id,
        "deleted_count": deleted_count,
        "created_count": len(created),
    }
