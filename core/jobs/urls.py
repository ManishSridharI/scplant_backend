from django.urls import re_path

from .api.TestAPI import Test
from .api.JobInferenceAPI import JobInference


app_name = 'jobs'


urlpatterns = [
    re_path('api/job_test/', Test, name='api_job_test'),
    re_path('api/job_inference/', JobInference, name='api_job_inference'),
]
