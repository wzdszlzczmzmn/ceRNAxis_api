from django.urls import path, include

from rest_framework.routers import DefaultRouter

from database.views import cerna_axis_views, cerna_axis_network_views

urlpatterns = [
    path('ceRNAAxis_table_filter_options/', cerna_axis_views.FilterOptionsView.as_view(),
         name='ceRNAAxis-table-filter-options'),
    path('ceRNAAxis_table_records/', cerna_axis_views.RNAInteractionSearchView.as_view(),
         name='ceRNAAxis-table-records'),
    path('ceRNAAxis_network_query/', cerna_axis_network_views.RNAInteractionNetworkView.as_view(),
         name='ceRNAAxis-network-query'),
]
