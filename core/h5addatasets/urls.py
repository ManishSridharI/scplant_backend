from django.urls import re_path

from .api.H5adDatasetGenerateAPI import H5adDatasetGenerate
from .api.H5adDatasetUploadAPI import H5adDatasetUpload
from .api.H5adDatasetQueryAPI import H5adDatasetQuery, H5adDatasetQueryPublic, H5adDatasetQueryUploadedAndPublic
from .api.H5adDatasetDeleteAPI import H5adDatasetDelete


app_name = 'h5addatasets'


urlpatterns = [
    re_path('api/h5ad_dataset_generate/', H5adDatasetGenerate, name='api_h5ad_dataset_generate'),
    re_path('api/h5ad_dataset_upload/', H5adDatasetUpload, name='api_h5ad_dataset_upload'),
    re_path('api/h5ad_dataset_query/', H5adDatasetQuery, name='api_h5ad_dataset_query'),
    re_path('api/h5ad_dataset_query_public/', H5adDatasetQueryPublic, name='api_h5ad_dataset_query_public'),
    re_path('api/h5ad_dataset_query_uploaded_and_public/', H5adDatasetQueryUploadedAndPublic, name='api_h5ad_dataset_query_uploaded_and_public'),
    re_path('api/h5ad_dataset_delete/', H5adDatasetDelete, name='api_h5ad_dataset_delete')
]
