from django.urls import path

from analysis.views import workflow_views, workflow_submit_views

urlpatterns = [
    path('immune_annotations/', workflow_views.WorkflowImmuneAnnotationListView.as_view(),
         name='workflow_immune_annotations'),
    path('immune_annotation/', workflow_views.WorkflowImmuneAnnotationDetailView.as_view(),
         name='workflow_immune_annotation_detail'),
    path('immune_annotation_download/', workflow_views.WorkflowImmuneAnnotationDownloadView.as_view(),
         name='workflow_immune_annotation_download'),
    path('custom_list_query_task_submit/', workflow_submit_views.CustomListQueryTaskSubmitView.as_view(),
         name='custom_list_query_task_submit'),
]
