from django.urls import re_path

from .api.TestAPI import Test
from .api.JobInferenceAPI import JobInference
from .api.JobInferenceQueryAPI import JobInferenceQuery
from .api.JobInferenceFileOutputQueryAPI import JobInferenceFileOutputQuery, JobInferenceFileOutputQueryByID
from .api.JobAnnotateAndPlotAPI import JobAnnotateAndPlot
from .api.JobAnnotateAndPlotQueryAPI import JobAnnotateAndPlotQuery
from .api.JobAnnotateAndPlotFileOutputQueryAPI import JobAnnotateAndPlotFileOutputQuery, JobAnnotateAndPlotFileOutputQueryByID


app_name = 'jobs'


urlpatterns = [
    re_path('api/job_test/', Test, name='api_job_test'),
    re_path('api/job_inference/', JobInference, name='api_job_inference'),
    re_path('api/job_inference_query/', JobInferenceQuery, name='api_job_inference_query'),
    re_path('api/job_inference_file_output_query/', JobInferenceFileOutputQuery, name='api_job_inference_file_output_query'),
    re_path('api/job_inference_file_output_query_by_id/', JobInferenceFileOutputQueryByID, name='api_job_inference_file_output_query_by_id'),
    re_path('api/job_annotate_and_plot/', JobAnnotateAndPlot, name='api_job_annotate_and_plot'),
    re_path('api/job_annotate_and_plot_query/', JobAnnotateAndPlotQuery, name='api_job_annotate_and_plot_query'),
    re_path('api/job_annotate_and_plot_file_output_query/', JobAnnotateAndPlotFileOutputQuery, name='api_job_annotate_and_plot_file_output_query'),
    re_path('api/job_annotate_and_plot_file_output_query_by_id/', JobAnnotateAndPlotFileOutputQueryByID, name='api_job_annotate_and_plot_file_output_query_by_id'),
]
