from django.urls import path

from database.views import cerna_axis_views, cerna_axis_network_views, dataset_views
from database.views.dataset_annotation_views import dataset_annotation_views, dataset_annotation_network_views, \
    dataset_annotation_axis_final_views, dataset_annotation_cmap_views, dataset_annotation_deg_volcano_views, \
    dataset_annotation_log2fc_correlation_views, dataset_annotation_survival_views, \
    dataset_annotation_deg_pathway_views, dataset_annotation_exp_correlation_views

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

    # Dataset File Download Views
    path('dataset_data_download/', dataset_views.DatasetDownloadView.as_view(),
         name='dataset-data-download'),
    path('tcga_annotation_download/', dataset_views.DatasetAnnotationDownloadView.as_view(),
         name='dataset-annotation-data-download'),
    path(
        "timedb_annotation_download/",
        dataset_views.TIMEDBAnnotationDownloadView.as_view(),
        name="timedb-annotation-download",
    ),

    # Sample Meta Views
    path('dataset_metadata/<str:dataset>/sample_meta/', dataset_views.DatasetSampleMetaView.as_view(),
         name='dataset-sample-meta'),
    path(
        "dataset_metadata/<str:dataset>/large_meta/", dataset_views.DatasetLargeMetaView.as_view(),
        name="dataset-large-meta",
    ),

    # Expression Data Views
    path('dataset_expression_genes/', dataset_views.DatasetExpressionGeneListView.as_view(),
         name='dataset-expression-genes'),
    path('dataset_expression/', dataset_views.DatasetExpressionDataView.as_view(),
         name='dataset-expression-data'),
    path(
        "dataset_expression/large_data/",
        dataset_views.DatasetLargeExpressionDataView.as_view(),
        name="dataset-large-expression-data",
    ),

    # DEG Volcano Views
    path('dataset_deg_volcano/', dataset_views.DatasetDEGVolcanoView.as_view(),
         name='dataset-deg-volcano'),
    path(
        "<str:dataset>/tisch2_deg_cluster_plot/",
        dataset_views.TISCH2DEGClusterPlotView.as_view(),
        name="tisch2-deg-cluster-plot",
    ),

    # Aliquot Files Download Views
    path('aliquot_expression_download_files/', dataset_views.DatasetAliquotExpressionFileListView.as_view(),
         name='aliquot-expression-download-files'),
    path('aliquot_expression_file_download/', dataset_views.DatasetAliquotExpressionFileDownloadView.as_view(),
         name='aliquot-expression-download-file'),

    # Annotations Available Views
    path('tcga_dataset_annotation_available/',
         dataset_annotation_views.TCGADatasetAnnotationAvailabilityView.as_view(),
         name='tcga-dataset-annotation-available'),
    path('timedb_dataset_annotation_available/',
         dataset_annotation_views.TIMEDBDatasetAnnotationAvailabilityView.as_view(),
         name='timedb-dataset-annotation-available'),

    # TIMEDB Annotation Options Views
    path("timedb_dataset_group_by_options/",
         dataset_annotation_views.TIMEDBAnnotationGroupByOptionsView.as_view(),
         name="timedb-annotation-group-by-options"),

    # Network Views
    path(
        "tcga_dataset_annotation_network/",
        dataset_annotation_network_views.TCGADatasetAnnotationNetworkView.as_view(),
        name="tcga_dataset_annotation_network",
    ),
    path(
        "timedb_dataset_annotation_network/",
        dataset_annotation_network_views.TIMEDBDatasetAnnotationNetworkView.as_view(),
        name="timedb_dataset_annotation_network",
    ),

    # Axis Final Views
    path(
        "tcga_dataset_annotation_axis_final/",
        dataset_annotation_axis_final_views.TCGADatasetAnnotationAxisFinalDataView.as_view(),
        name="tcga_dataset_annotation_axis_final",
    ),
    path(
        "timedb_dataset_annotation_axis_final/",
        dataset_annotation_axis_final_views.TIMEDBDatasetAnnotationAxisFinalDataView.as_view(),
        name="timedb_dataset_annotation_axis_final",
    ),

    # CMAP Views
    path(
        "tcga_dataset_annotation_cmap/",
        dataset_annotation_cmap_views.TCGADatasetAnnotationCMapResultView.as_view(),
        name="tcga_dataset_annotation_cmap",
    ),
    path(
        "timedb_dataset_annotation_cmap/",
        dataset_annotation_cmap_views.TIMEDBDatasetAnnotationCMapResultView.as_view(),
        name="timedb_dataset_annotation_cmap",
    ),

    # DEG Volcano
    path(
        "tcga_dataset_annotation_deg_volcano/",
        dataset_annotation_deg_volcano_views.TCGADatasetAnnotationDEGVolcanoView.as_view(),
        name="tcga_dataset_annotation_deg_volcano",
    ),
    path(
        "timedb_dataset_annotation_deg_volcano/",
        dataset_annotation_deg_volcano_views.TIMEDBDatasetAnnotationDEGVolcanoView.as_view(),
        name="timedb_dataset_annotation_deg_volcano",
    ),

    # Log2FC Correlation Views
    path(
        "tcga_dataset_annotation_log2fc_correlation/",
        dataset_annotation_log2fc_correlation_views.TCGADatasetAnnotationLog2FCCorrelationView.as_view(),
        name="tcga_dataset_annotation_log2fc_correlation",
    ),
    path(
        "timedb_dataset_annotation_log2fc_correlation/",
        dataset_annotation_log2fc_correlation_views.TIMEDBDatasetAnnotationLog2FCCorrelationView.as_view(),
        name="timedb_dataset_annotation_log2fc_correlation",
    ),

    # Exp Correlation Views
    path(
        "tcga_dataset_annotation_exp_correlation_options/",
        dataset_annotation_exp_correlation_views.TCGADatasetAnnotationExpCorrelationOptionsView.as_view(),
        name="tcga_dataset_annotation_exp_correlation_options",
    ),
    path(
        "tcga_dataset_annotation_exp_correlation_plot_data/",
        dataset_annotation_exp_correlation_views.TCGADatasetAnnotationExpCorrelationPlotDataView.as_view(),
        name="tcga_dataset_annotation_exp_correlation_plot_data",
    ),
    path(
        "timedb_dataset_annotation_exp_correlation_options/",
        dataset_annotation_exp_correlation_views.TIMEDBDatasetAnnotationExpCorrelationOptionsView.as_view(),
        name="timedb_dataset_annotation_exp_correlation_options",
    ),
    path(
        "timedb_dataset_annotation_exp_correlation_plot_data/",
        dataset_annotation_exp_correlation_views.TIMEDBDatasetAnnotationExpCorrelationPlotDataView.as_view(),
        name="timedb_dataset_annotation_exp_correlation_plot_data",
    ),

    # Survival KM Views
    path(
        "tcga_dataset_annotation_survival_km/",
        dataset_annotation_survival_views.TCGADatasetAnnotationSurvivalKMDataView.as_view(),
        name="tcga_dataset_annotation_survival_km",
    ),
    path(
        "timedb_dataset_annotation_survival_km/",
        dataset_annotation_survival_views.TIMEDBDatasetAnnotationSurvivalKMDataView.as_view(),
        name="timedb_dataset_annotation_survival_km",
    ),

    # DEG Pathway Views
    path(
        "tcga_dataset_annotation_deg_pathway/",
        dataset_annotation_deg_pathway_views.TCGADatasetAnnotationDEGPathwayView.as_view(),
        name="tcga_dataset_annotation_deg_pathway",
    ),
    path(
        "timedb_dataset_annotation_deg_pathway/",
        dataset_annotation_deg_pathway_views.TIMEDBDatasetAnnotationDEGPathwayView.as_view(),
        name="timedb_dataset_annotation_deg_pathway",
    ),
]
