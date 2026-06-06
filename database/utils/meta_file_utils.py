from django.conf import settings


def get_dataset_meta_file(dataset: str):
    dataset_prefix = "_".join(dataset.split("_")[:2])
    return settings.DATASET_BASE_DIR / f"{dataset_prefix}_meta.csv"
