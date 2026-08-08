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
)
from database.utils.dataset_annotation_utils.metadata_utils import (
    build_dataset_annotation_metadata,
    build_tcga_dataset_annotation_availability,
    build_timedb_dataset_annotation_availability,
    DEFAULT_DEG_METHOD,
)
from database.utils.dataset_annotation_utils.scst_metadata_utils import build_scst_dataset_annotation_availability
from database.utils.dataset_annotation_utils.scst_path_utils import get_scst_data_type_query


class TCGADatasetAnnotationAvailabilityView(APIView):
    source = "TCGA"
    annotation_root_setting_name = "TCGA_DATASET_ANNOTATIONS_DIR"
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

    def get_deg_method(self, request) -> str:
        deg_method = request.query_params.get(
            "deg_method",
            self.default_deg_method,
        )

        return str(
            deg_method or self.default_deg_method
        ).strip() or self.default_deg_method

    def get(self, request):
        try:
            dataset_name = get_dataset_query_name(request)
            deg_method = self.get_deg_method(request)

            annotation_dir_name = resolve_tcga_annotation_dir_name(
                dataset_name
            )

            annotation_dir = resolve_dataset_annotation_dir(
                annotation_root_dir=self.get_annotation_root_dir(),
                annotation_dir_name=annotation_dir_name,
            )

            availability = build_tcga_dataset_annotation_availability(
                annotation_dir=annotation_dir,
                file_prefix=annotation_dir_name,
                deg_method=deg_method,
            )

            return Response(
                {
                    "success": True,
                    "source": self.source,
                    "dataset_name": dataset_name,
                    **availability,
                },
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


class TIMEDBDatasetAnnotationAvailabilityView(APIView):
    source = "TIMEDB"
    annotation_root_setting_name = "TIMEDB_DATASET_ANNOTATIONS_DIR"
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

    def get_deg_method(self, request) -> str:
        deg_method = request.query_params.get(
            "deg_method",
            self.default_deg_method,
        )

        return str(
            deg_method or self.default_deg_method
        ).strip() or self.default_deg_method

    def get(self, request):
        try:
            dataset_name = get_dataset_query_name(request)
            deg_method = self.get_deg_method(request)

            availability = build_timedb_dataset_annotation_availability(
                annotation_root_dir=self.get_annotation_root_dir(),
                dataset_name=dataset_name,
                deg_method=deg_method,
            )

            return Response(
                {
                    "success": True,
                    "source": self.source,
                    "dataset_name": dataset_name,
                    **availability,
                },
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


class SCSTDatasetAnnotationAvailabilityView(APIView):
    """
    Return SC/ST Dataset Annotation availability.

    Query parameters:
        dataset:
            Dataset name, e.g. ALL_GSE154109.

        data_type:
            sc | st

    Current scope:
        - configured group-by options
        - group-by result-directory availability
        - unique group values and counts from dataset meta CSV

    Per-group-value visualization availability is deferred.
    """

    source = "SCST"

    def get(self, request):
        try:
            dataset_name = get_dataset_query_name(
                request
            )

            data_type = get_scst_data_type_query(
                request
            )

            availability = (
                build_scst_dataset_annotation_availability(
                    dataset_name=dataset_name,
                    data_type=data_type,
                )
            )

            return Response(
                {
                    "success": True,
                    "source": self.source,
                    "data_type": data_type,
                    "dataset_name": dataset_name,
                    **availability,
                },
                status=status.HTTP_200_OK,
            )

        except DatasetAnnotationInputError as exc:
            return Response(
                {
                    "success": False,
                    "detail": str(exc),
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        except FileNotFoundError as exc:
            return Response(
                {
                    "success": False,
                    "detail": str(exc),
                },
                status=status.HTTP_404_NOT_FOUND,
            )

        except DatasetAnnotationPathError as exc:
            return Response(
                {
                    "success": False,
                    "detail": str(exc),
                },
                status=(
                    status.HTTP_500_INTERNAL_SERVER_ERROR
                ),
            )

        except Exception as exc:
            print(traceback.format_exc())

            return Response(
                {
                    "success": False,
                    "detail": (
                        f"Server error: {str(exc)}"
                    ),
                },
                status=(
                    status.HTTP_500_INTERNAL_SERVER_ERROR
                ),
            )

