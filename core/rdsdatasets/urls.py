from django.urls import re_path

from .api.RdsDatasetGenerateAPI import RdsDatasetGenerate
from .api.RdsDatasetUploadAPI import RdsDatasetUpload
from .api.RdsDatasetQueryAPI import RdsDatasetQuery, RdsDatasetQueryPublic, RdsDatasetQueryUploadedAndPublic
from .api.RdsDatasetDeleteAPI import RdsDatasetDelete


app_name = 'rdsdatasets'


urlpatterns = [
    re_path('api/rds_dataset_generate/', RdsDatasetGenerate, name='api_rds_dataset_generate'),
    re_path('api/rds_dataset_upload/', RdsDatasetUpload, name='api_rds_dataset_upload'),
    re_path('api/rds_dataset_query/', RdsDatasetQuery, name='api_rds_dataset_query'),
    re_path('api/rds_dataset_query_public/', RdsDatasetQueryPublic, name='api_rds_dataset_query_public'),
    re_path('api/rds_dataset_query_uploaded_and_public/', RdsDatasetQueryUploadedAndPublic, name='api_rds_dataset_query_uploaded_and_public'),
    re_path('api/rds_dataset_delete/', RdsDatasetDelete, name='api_rds_dataset_delete')
]
