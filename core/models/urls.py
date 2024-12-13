from django.urls import re_path

from .api.ModelUploadAPI import ModelUpload
from .api.ModelQueryAPI import ModelQuery, ModelQueryPublic, ModelQueryUploadedAndPublic
from .api.ModelDeleteAPI import ModelDelete


app_name = 'models'


urlpatterns = [
    re_path('api/model_upload/', ModelUpload, name='api_model_upload'),
    re_path('api/model_query/', ModelQuery, name='api_model_query'),
    re_path('api/model_query_public/', ModelQueryPublic, name='api_model_query_public'),
    re_path('api/model_query_uploaded_and_public/', ModelQueryUploadedAndPublic, name='api_model_query_uploaded_and_public'),
    re_path('api/model_delete/', ModelDelete, name='api_model_delete')
]
