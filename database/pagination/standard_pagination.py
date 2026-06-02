from rest_framework.pagination import PageNumberPagination
from rest_framework.response import Response


class StandardPageNumberPagination(PageNumberPagination):
    page_query_param = "page"
    page_size_query_param = "page_size"
    page_size = 10
    max_page_size = 100

    def get_page_number(self, request, paginator):
        if request.method == "POST":
            return request.data.get(self.page_query_param, 1)

        return request.query_params.get(self.page_query_param, 1)

    def get_page_size(self, request):
        if request.method == "POST":
            page_size = request.data.get(self.page_size_query_param)

            if page_size is not None:
                try:
                    page_size = int(page_size)
                except (TypeError, ValueError):
                    return self.page_size

                return min(page_size, self.max_page_size)

        return super().get_page_size(request)

    def get_paginated_response(self, data):
        return Response({
            "count": self.page.paginator.count,
            "page": self.page.number,
            "page_size": self.get_page_size(self.request),
            "num_pages": self.page.paginator.num_pages,
            "results": data,
        })
