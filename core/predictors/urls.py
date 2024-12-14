from django.urls import re_path

from .api.PredictorUploadAPI import PredictorUpload
from .api.PredictorQueryAPI import PredictorQuery, PredictorQueryPublic, PredictorQueryUploadedAndPublic
from .api.PredictorDeleteAPI import PredictorDelete


app_name = 'predictors'


urlpatterns = [
    re_path('api/predictor_upload/', PredictorUpload, name='api_predictor_upload'),
    re_path('api/predictor_query/', PredictorQuery, name='api_predictor_query'),
    re_path('api/predictor_query_public/', PredictorQueryPublic, name='api_predictor_query_public'),
    re_path('api/predictor_query_uploaded_and_public/', PredictorQueryUploadedAndPublic, name='api_predictor_query_uploaded_and_public'),
    re_path('api/predictor_delete/', PredictorDelete, name='api_predictor_delete')
]
