import pandas as pd
from django.http import FileResponse

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from analysis.utils.immune_annotation_path_utils import get_immune_annotation_raw_dir, IMMUNE_FILE_PREFIX, \
    IMMUNE_FILE_SUFFIX, get_immune_annotation_file_path, get_label_from_map_info, ImmuneAnnotationPathError, \
    validate_immune_annotation_file


class WorkflowImmuneAnnotationListView(APIView):
    def get(self, request):
        raw_dir = get_immune_annotation_raw_dir()

        if not raw_dir.exists() or not raw_dir.is_dir():
            return Response(
                {"detail": "Immune annotation raw directory not found."},
                status=status.HTTP_404_NOT_FOUND,
            )

        results = []

        for file_path in sorted(raw_dir.glob(f"{IMMUNE_FILE_PREFIX}*{IMMUNE_FILE_SUFFIX}")):
            if not file_path.is_file():
                continue

            stem = file_path.stem
            label = stem.removeprefix(IMMUNE_FILE_PREFIX)

            if not label:
                continue

            results.append(
                {
                    "label": label,
                    "value": stem,
                }
            )

        return Response(
            {
                "count": len(results),
                "results": results,
            },
            status=status.HTTP_200_OK,
        )


class WorkflowImmuneAnnotationDetailView(APIView):
    def get(self, request):
        map_info = request.query_params.get("map_info")

        try:
            file_path = get_immune_annotation_file_path(map_info)
            label = get_label_from_map_info(map_info)
        except ImmuneAnnotationPathError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )

        if not file_path.exists() or not file_path.is_file():
            return Response(
                {
                    "detail": "Immune annotation file not found.",
                    "map_info": map_info,
                },
                status=status.HTTP_404_NOT_FOUND,
            )

        try:
            df = pd.read_csv(file_path)
        except Exception as e:
            return Response(
                {"detail": f"Failed to read immune annotation file: {e}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        df = df.fillna("")

        return Response(
            {
                "map_info": map_info,
                "label": label,
                "count": int(df.shape[0]),
                "columns": df.columns.tolist(),
                "results": df.to_dict("records"),
            },
            status=status.HTTP_200_OK,
        )


class WorkflowImmuneAnnotationDownloadView(APIView):
    def get(self, request):
        map_info = request.query_params.get("map_info")

        try:
            file_path = validate_immune_annotation_file(map_info)
        except ImmuneAnnotationPathError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )
        except FileNotFoundError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_404_NOT_FOUND,
            )

        download_filename = f"{map_info}{IMMUNE_FILE_SUFFIX}"

        return FileResponse(
            open(file_path, "rb"),
            as_attachment=True,
            filename=download_filename,
            content_type="text/csv",
        )
