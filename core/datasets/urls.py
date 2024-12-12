from django.urls import re_path

from .api.DatasetUploadAPI import DatasetUpload
from .api.DatasetQueryAPI import DatasetQuery
from .api.DatasetDeleteAPI import DatasetDelete


app_name = 'datasets'


urlpatterns = [
    re_path('api/dataset_upload/', DatasetUpload, name='api_dataset_upload'),
    re_path('api/dataset_query/', DatasetQuery, name='api_dataset_query'),
    re_path('api/dataset_delete/', DatasetDelete, name='api_dataset_delete')
]
