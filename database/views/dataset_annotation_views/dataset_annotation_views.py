import traceback

from django.conf import settings
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from database.utils.dataset_annotation_utils.path_utils import (
    DatasetAnnotationInputError,
    DatasetAnnotationPathError,
    get_dataset_query_name,
    resolve_dataset_annotation_dir,
    resolve_tcga_annotation_dir_name,
    resolve_timedb_annotation_dir_name,
)
from database.utils.dataset_annotation_utils.metadata_utils import (
    build_dataset_annotation_metadata,
    DEFAULT_DEG_METHOD,
)


class BaseDatasetAnnotationAvailabilityView(APIView):
    source = None
    annotation_root_setting_name = None
    annotation_dir_name_resolver = None

    default_deg_method = DEFAULT_DEG_METHOD

    def get_annotation_root_dir(self):
        annotation_root_dir = getattr(
            settings,
            self.annotation_root_setting_name,
            None,
        )

        if not annotation_root_dir:
            raise DatasetAnnotationPathError(
                f"{self.annotation_root_setting_name} is not configured."
            )

        return annotation_root_dir

    def get_annotation_file_prefix(
        self,
        *,
        dataset_name: str,
        annotation_dir_name: str,
    ) -> str:
        return annotation_dir_name

    def get(self, request):
        try:
            dataset_name = get_dataset_query_name(request)

            deg_method = str(
                request.query_params.get(
                    "deg_method",
                    self.default_deg_method,
                )
            ).strip()

            annotation_dir_name = self.annotation_dir_name_resolver(
                dataset_name
            )

            annotation_dir = resolve_dataset_annotation_dir(
                annotation_root_dir=self.get_annotation_root_dir(),
                annotation_dir_name=annotation_dir_name,
            )

            file_prefix = self.get_annotation_file_prefix(
                dataset_name=dataset_name,
                annotation_dir_name=annotation_dir_name,
            )

            available = annotation_dir.exists() and annotation_dir.is_dir()

            response_data = {
                "success": True,
                "source": self.source,
                "dataset_name": dataset_name,
                "annotation_dir_name": annotation_dir_name,
                "annotation_file_prefix": file_prefix,
                "available": available,
            }

            if available:
                response_data.update(
                    build_dataset_annotation_metadata(
                        annotation_dir=annotation_dir,
                        file_prefix=file_prefix,
                        source=self.source,
                        deg_method=deg_method,
                    )
                )
            else:
                response_data.update(
                    {
                        "network_source_task_type": None,
                        "deg_method": deg_method,
                        "use_padj": True,
                        "cutoffs": {},
                        "available_deg_rna_types": [],
                        "available_deg_scopes": [],
                        "available_background_types": [],
                    }
                )

            return Response(
                response_data,
                status=status.HTTP_200_OK,
            )

        except DatasetAnnotationInputError as e:
            return Response(
                {
                    "success": False,
                    "detail": str(e),
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        except DatasetAnnotationPathError as e:
            return Response(
                {
                    "success": False,
                    "detail": str(e),
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        except Exception as e:
            print(traceback.format_exc())

            return Response(
                {
                    "success": False,
                    "detail": f"Server error: {str(e)}",
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )


class TCGADatasetAnnotationAvailabilityView(
    BaseDatasetAnnotationAvailabilityView
):
    source = "TCGA"
    annotation_root_setting_name = "TCGA_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_tcga_annotation_dir_name
    )


class TIMEDBDatasetAnnotationAvailabilityView(
    BaseDatasetAnnotationAvailabilityView
):
    source = "TIMEDB"
    annotation_root_setting_name = "TIMEDB_DATASET_ANNOTATIONS_DIR"
    annotation_dir_name_resolver = staticmethod(
        resolve_timedb_annotation_dir_name
    )
