from django.urls import re_path

from .api.TenxfbcmDatasetGenerateAPI import TenxfbcmDatasetGenerate
from .api.TenxfbcmDatasetUploadAPI import TenxfbcmDatasetUpload
from .api.TenxfbcmDatasetQueryAPI import TenxfbcmDatasetQuery, TenxfbcmDatasetQueryPublic, TenxfbcmDatasetQueryUploadedAndPublic
from .api.TenxfbcmDatasetDeleteAPI import TenxfbcmDatasetDelete


app_name = 'tenxfeaturebcmatrixdatasets'


urlpatterns = [
    re_path('api/tenxfbcm_dataset_generate/', TenxfbcmDatasetGenerate, name='api_tenxfbcm_dataset_generate'),
    re_path('api/tenxfbcm_dataset_upload/', TenxfbcmDatasetUpload, name='api_tenxfbcm_dataset_upload'),
    re_path('api/tenxfbcm_dataset_query/', TenxfbcmDatasetQuery, name='api_tenxfbcm_dataset_query'),
    re_path('api/tenxfbcm_dataset_query_public/', TenxfbcmDatasetQueryPublic, name='api_tenxfbcm_dataset_query_public'),
    re_path('api/tenxfbcm_dataset_query_uploaded_and_public/', TenxfbcmDatasetQueryUploadedAndPublic, name='api_tenxfbcm_dataset_query_uploaded_and_public'),
    re_path('api/tenxfbcm_dataset_delete/', TenxfbcmDatasetDelete, name='api_tenxfbcm_dataset_delete')
]
