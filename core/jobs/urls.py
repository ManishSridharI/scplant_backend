from django.urls import re_path

from .api.TestAPI import Test
from .api.JobPredictorInferenceAPI import JobPredictorInference


app_name = 'jobs'


urlpatterns = [
    re_path('api/test/', Test, name='api_test'),
    re_path('api/job_predictor_inference/', JobPredictorInference, name='api_job_predictor_inference'),
]
