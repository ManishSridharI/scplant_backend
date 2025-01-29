from django.urls import re_path

from .api.PredDatasetGenerateAPI import PredDatasetGenerate
from .api.PredDatasetUploadAPI import PredDatasetUpload
from .api.PredDatasetQueryAPI import PredDatasetQuery, PredDatasetQueryPublic, PredDatasetQueryUploadedAndPublic
from .api.PredDatasetDeleteAPI import PredDatasetDelete


app_name = 'preddatasets'


urlpatterns = [
    re_path('api/pred_dataset_generate/', PredDatasetGenerate, name='api_pred_dataset_generate'),
    re_path('api/pred_dataset_upload/', PredDatasetUpload, name='api_pred_dataset_upload'),
    re_path('api/pred_dataset_query/', PredDatasetQuery, name='api_pred_dataset_query'),
    re_path('api/pred_dataset_query_public/', PredDatasetQueryPublic, name='api_pred_dataset_query_public'),
    re_path('api/pred_dataset_query_uploaded_and_public/', PredDatasetQueryUploadedAndPublic, name='api_pred_dataset_query_uploaded_and_public'),
    re_path('api/pred_dataset_delete/', PredDatasetDelete, name='api_pred_dataset_delete')
]
