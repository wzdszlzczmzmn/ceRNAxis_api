from django.urls import path

from analysis.views import workflow_views, workflow_submit_views, workflow_query_views, workflow_detail_views, \
    workflow_demo_view
from analysis.views.workflow_detail_each_views import workflow_uploaded_file_views, workflow_network_views, \
    workflow_axis_final_views, workflow_cmap_views, workflow_log2fc_correlation_views, workflow_deg_volcano_views, \
    workflow_exp_correlation_views, workflow_survival_views, workflow_deg_pathway_views

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

    # Task Network Views
    path('custom_list_query_network/', workflow_network_views.CustomListQueryTaskNetworkView.as_view(),
         name='custom_list_query_network'),
    path('paired_cohort_task_network/', workflow_network_views.PairedCohortTaskNetworkView.as_view(),
         name='paired_cohort_task_network'),
    path('hybrid_reference_task_network/', workflow_network_views.HybridReferenceTaskNetworkView.as_view(),
         name='hybrid_reference_task_network'),

    path('paired_cohort_task_submit/', workflow_submit_views.PairedCohortTaskSubmitView.as_view(),
         name='paired_cohort_task_submit'),

    # Uploaded File Download Views
    path('paired_cohort_uploaded_file_download/',
         workflow_uploaded_file_views.PairedCohortUploadedFileDownloadView.as_view(),
         name='paired_cohort_uploaded_file_download'),
    path('hybrid_reference_uploaded_file_download/',
         workflow_uploaded_file_views.HybridReferenceUploadedFileDownloadView.as_view(),
         name='hybrid_reference_uploaded_file_download'
         ),

    # DEG Volcano Views
    path('paired_cohort_deg_volcano/', workflow_deg_volcano_views.PairedCohortDEGVolcanoView.as_view(),
         name='paired_cohort_deg_volcano'),
    path('hybrid_reference_deg_volcano/', workflow_deg_volcano_views.HybridReferenceDEGVolcanoView.as_view(),
         name='hybrid_reference_deg_volcano'),

    # Log2FC Correlation Views
    path('paired_cohort_correlation/', workflow_log2fc_correlation_views.PairedCohortLog2FCCorrelationView.as_view(),
         name='paired_cohort_correlation'),
    path('hybrid_reference_correlation/',
         workflow_log2fc_correlation_views.HybridReferenceLog2FCCorrelationView.as_view(),
         name='hybrid_reference_correlation'),

    # Exp Correlation Views
    path('paired_cohort_exp_correlation_options/',
         workflow_exp_correlation_views.PairedCohortExpCorrelationOptionsView.as_view(),
         name='paired_cohort_exp_correlation_options'),
    path('paired_cohort_exp_correlation_plot_data/',
         workflow_exp_correlation_views.PairedCohortExpCorrelationPlotDataView.as_view(),
         name='paired_cohort_exp_correlation_plot_data'),
    path('hybrid_reference_exp_correlation_options/',
         workflow_exp_correlation_views.HybridReferenceExpCorrelationOptionsView.as_view(),
         name='hybrid_reference_exp_correlation_options'),
    path('hybrid_reference_exp_correlation_plot_data/',
         workflow_exp_correlation_views.HybridReferenceExpCorrelationPlotDataView.as_view(),
         name='hybrid_reference_exp_correlation_plot_data'),

    # Task Demo
    path('paired_cohort_demo_info/', workflow_demo_view.PairedCohortDemoInfoView.as_view(),
         name='paired_cohort_demo_info'),
    path('paired_cohort_demo_sample_meta/', workflow_demo_view.PairedCohortDemoSampleMetaView.as_view(),
         name='paired_cohort_demo_sample_meta'),
    path('paired_cohort_demo_expression_gene_list/',
         workflow_demo_view.PairedCohortDemoExpressionGeneListView.as_view(),
         name='paired_cohort_demo_expression_gene_list'),
    path('paired_cohort_demo_expression_data/', workflow_demo_view.PairedCohortDemoExpressionDataView.as_view(),
         name='paired_cohort_demo_expression_data'),
    path('hybrid_reference_demo_info/', workflow_demo_view.HybridReferenceDemoInfoView.as_view(),
         name='hybrid_reference_demo_info'),
    path('hybrid_reference_demo_sample_meta/', workflow_demo_view.HybridReferenceDemoSampleMetaView.as_view(),
         name='hybrid_reference_demo_sample_meta'),
    path('hybrid_reference_demo_expression_gene_list/',
         workflow_demo_view.HybridReferenceDemoExpressionGeneListView.as_view(),
         name='hybrid_reference_demo_expression_gene_list'),
    path('hybrid_reference_demo_expression_data/', workflow_demo_view.HybridReferenceDemoExpressionDataView.as_view(),
         name='hybrid_reference_demo_expression_data'),

    # Demo Input Download Views
    path('paired_cohort_demo_download_data/', workflow_demo_view.PairedCohortDemoDataDownloadView.as_view(),
         name='paired_cohort_demo_download_data'),
    path('hybrid_reference_demo_files_download/', workflow_demo_view.HybridReferenceDemoDataDownloadView.as_view(),
         name='hybrid_reference_demo_files_download'),

    path('workflow_task_result_download/', workflow_detail_views.WorkflowTaskResultDownloadView.as_view(),
         name='workflow_task_result_download'),

    # Run Demo Views
    path('custom_list_query_run_demo/', workflow_demo_view.CustomListQueryDemoRunView.as_view(),
         name='custom_list_query_run_demo'),
    path('paired_cohort_run_demo/', workflow_demo_view.PairedCohortDemoRunView.as_view(),
         name='paired_cohort_run_demo'),
    path("hybrid_reference_run_demo/", workflow_demo_view.HybridReferenceDemoRunView.as_view(),
         name="hybrid_reference_demo_run"),

    # Axis Final Views
    path('paired_cohort_axis_final/', workflow_axis_final_views.PairedCohortAxisFinalDataView.as_view(),
         name='paired_cohort_axis_final'),
    path('hybrid_reference_axis_final/', workflow_axis_final_views.HybridReferenceAxisFinalDataView.as_view(),
         name='hybrid_reference_axis_final'),

    # Survival Views
    path('paired_cohort_survival_km/', workflow_survival_views.PairedCohortSurvivalKMDataView.as_view(),
         name='paired_cohort_survival_km'),
    path('hybrid_reference_survival_km/', workflow_survival_views.HybridReferenceSurvivalKMDataView.as_view(),
         name='hybrid_reference_survival_km'),

    path('hybrid_reference_task_submit/', workflow_submit_views.HybridReferenceTaskSubmitView.as_view(),
         name='hybrid_reference_task_submit'),

    # DEG Pathway Views
    path('paired_cohort_deg_pathway/', workflow_deg_pathway_views.PairedCohortDEGPathwayView.as_view(),
         name='paired_cohort_deg_pathway'),
    path('hybrid_reference_deg_pathway/', workflow_deg_pathway_views.HybridReferenceDEGPathwayView.as_view(),
         name='hybrid_reference_deg_pathway'),

    # CMap Views
    path('paired_cohort_cmap_result/', workflow_cmap_views.PairedCohortCMapResultView.as_view(),
         name='paired_cohort_cmap_result'),
    path('hybrid_reference_cmap_result/', workflow_cmap_views.HybridReferenceCMapResultView.as_view(),
         name='hybrid_reference_cmap_result'),
]
