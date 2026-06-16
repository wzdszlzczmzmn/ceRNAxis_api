from django.urls import path

from analysis.views import workflow_views, workflow_submit_views, workflow_query_views, workflow_detail_views

urlpatterns = [
    path('immune_annotations/', workflow_views.WorkflowImmuneAnnotationListView.as_view(),
         name='workflow_immune_annotations'),
    path('immune_annotation/', workflow_views.WorkflowImmuneAnnotationDetailView.as_view(),
         name='workflow_immune_annotation_detail'),
    path('immune_annotation_download/', workflow_views.WorkflowImmuneAnnotationDownloadView.as_view(),
         name='workflow_immune_annotation_download'),
    path('custom_list_query_task_submit/', workflow_submit_views.CustomListQueryTaskSubmitView.as_view(),
         name='custom_list_query_task_submit'),
    path('query_task/', workflow_query_views.QueryTaskView.as_view(), name='query_task'),
    path('custom_list_query_network/', workflow_detail_views.WorkflowRNAInteractionNetworkView.as_view(),
         name='custom_list_query_network'),
    path('paired_cohort_task_submit/', workflow_submit_views.PairedCohortTaskSubmitView.as_view(),
         name='paired_cohort_task_submit'),
]
