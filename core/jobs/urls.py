from django.urls import re_path

from .api.TestAPI import Test
from .api.JobInferenceAPI import JobInference
from .api.JobInferenceQueryAPI import JobInferenceQuery
from .api.JobInferenceFileOutputQueryAPI import JobInferenceFileOutputQuery, JobInferenceFileOutputQueryByID
from .api.JobAnnotateAndPlotAPI import JobAnnotateAndPlot
from .api.JobAnnotateAndPlotQueryAPI import JobAnnotateAndPlotQuery
from .api.JobAnnotateAndPlotFileOutputQueryAPI import JobAnnotateAndPlotFileOutputQuery, JobAnnotateAndPlotFileOutputQueryByID
from .api.JobTreatmentVsControlAPI import JobTreatmentVsControl
from .api.JobTreatmentVsControlQueryAPI import JobTreatmentVsControlQuery
from .api.JobTreatmentVsControlFileOutputQueryAPI import JobTreatmentVsControlFileOutputQuery, JobTreatmentVsControlFileOutputQueryByID
from .api.JobCompareCellTypeDistAPI import JobCompareCellTypeDist
from .api.JobCompareCellTypeDistQueryAPI import JobCompareCellTypeDistQuery
from .api.JobCompareCellTypeDistFileOutputQueryAPI import JobCompareCellTypeDistFileOutputQuery, JobCompareCellTypeDistFileOutputQueryByID


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
    re_path('api/job_treatment_vs_control/', JobTreatmentVsControl, name='api_job_treatment_vs_control'),
    re_path('api/job_treatment_vs_control_query/', JobTreatmentVsControlQuery, name='api_job_treatment_vs_control_query'),
    re_path('api/job_treatment_vs_control_file_output_query/', JobTreatmentVsControlFileOutputQuery, name='api_job_treatment_vs_control_file_output_query'),
    re_path('api/job_treatment_vs_control_file_output_query_by_id/', JobTreatmentVsControlFileOutputQueryByID, name='api_job_treatment_vs_control_file_output_query_by_id'),
    re_path('api/job_compare_cell_type_dist/', JobCompareCellTypeDist, name='api_job_compare_cell_type_dist'),
    re_path('api/job_compare_cell_type_dist_query/', JobCompareCellTypeDistQuery, name='api_job_compare_cell_type_dist_query'),
    re_path('api/job_compare_cell_type_dist_file_output_query/', JobCompareCellTypeDistFileOutputQuery, name='api_job_compare_cell_type_dist_file_output_query'),
    re_path('api/job_compare_cell_type_dist_file_output_query_by_id/', JobCompareCellTypeDistFileOutputQueryByID, name='api_job_compare_cell_type_dist_file_output_query_by_id'),
]
