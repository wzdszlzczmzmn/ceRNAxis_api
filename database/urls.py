from django.urls import path

from database.views import cerna_axis_views, cerna_axis_network_views, dataset_views

urlpatterns = [
    path('ceRNAAxis_table_filter_options/', cerna_axis_views.FilterOptionsView.as_view(),
         name='ceRNAAxis-table-filter-options'),
    path('ceRNAAxis_table_records/', cerna_axis_views.RNAInteractionSearchView.as_view(),
         name='ceRNAAxis-table-records'),
    path('ceRNAAxis_network_query/', cerna_axis_network_views.RNAInteractionNetworkView.as_view(),
         name='ceRNAAxis-network-query'),
    path('dataset_metadata/', dataset_views.DatasetMetadataListView.as_view(), name="dataset-metadata-list"),
    path('dataset_metadata/<str:dataset>/', dataset_views.DatasetMetadataDetailView.as_view(),
         name='dataset-metadata-detail'),
    path('dataset_metadata/<str:dataset>/sample_meta/', dataset_views.DatasetSampleMetaView.as_view(),
         name='dataset-sample-meta'),
    path('dataset_expression_genes/', dataset_views.DatasetExpressionGeneListView.as_view(),
         name='dataset-expression-genes'),
    path('dataset_expression/', dataset_views.DatasetExpressionDataView.as_view(),
         name='dataset-expression-data'),
    path('dataset_deg_volcano/', dataset_views.DatasetDEGVolcanoView.as_view(),
         name='dataset-deg-volcano'),
    path('aliquot_expression_download_files/', dataset_views.DatasetAliquotExpressionFileListView.as_view(),
         name='aliquot-expression-download-files'),
    path('aliquot_expression_file_download/', dataset_views.DatasetAliquotExpressionFileDownloadView.as_view(),
         name='aliquot-expression-download-file'),
]
