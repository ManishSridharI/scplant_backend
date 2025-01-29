from django.urls import re_path

from .api.JobAnnotateAndPlotAPI import JobAnnotateAndPlot
from .api.JobAnnotateAndPlotQueryAPI import JobAnnotateAndPlotQuery
from .api.JobAnnotateAndPlotFileOutputQueryAPI import JobAnnotateAndPlotFileOutputQuery, JobAnnotateAndPlotFileOutputQueryByID
from .api.JobTreatmentVsControlAPI import JobTreatmentVsControl
from .api.JobTreatmentVsControlQueryAPI import JobTreatmentVsControlQuery
from .api.JobTreatmentVsControlFileOutputQueryAPI import JobTreatmentVsControlFileOutputQuery, JobTreatmentVsControlFileOutputQueryByID
from .api.JobConvertRdsToH5adAPI import JobConvertRdsToH5ad
from .api.JobConvertRdsToH5adQueryAPI import JobConvertRdsToH5adQuery
from .api.JobConvertRdsToH5adFileOutputQueryAPI import JobConvertRdsToH5adFileOutputQuery, JobConvertRdsToH5adFileOutputQueryByID


app_name = 'jobs'


urlpatterns = [
    re_path('api/job_annotate_and_plot/', JobAnnotateAndPlot, name='api_job_annotate_and_plot'),
    re_path('api/job_annotate_and_plot_query/', JobAnnotateAndPlotQuery, name='api_job_annotate_and_plot_query'),
    re_path('api/job_annotate_and_plot_file_output_query/', JobAnnotateAndPlotFileOutputQuery, name='api_job_annotate_and_plot_file_output_query'),
    re_path('api/job_annotate_and_plot_file_output_query_by_id/', JobAnnotateAndPlotFileOutputQueryByID, name='api_job_annotate_and_plot_file_output_query_by_id'),
    re_path('api/job_treatment_vs_control/', JobTreatmentVsControl, name='api_job_treatment_vs_control'),
    re_path('api/job_treatment_vs_control_query/', JobTreatmentVsControlQuery, name='api_job_treatment_vs_control_query'),
    re_path('api/job_treatment_vs_control_file_output_query/', JobTreatmentVsControlFileOutputQuery, name='api_job_treatment_vs_control_file_output_query'),
    re_path('api/job_treatment_vs_control_file_output_query_by_id/', JobTreatmentVsControlFileOutputQueryByID, name='api_job_treatment_vs_control_file_output_query_by_id'),
    re_path('api/job_convert_rds_to_h5ad/', JobConvertRdsToH5ad, name='api_job_convert_rds_to_h5ad'),
    re_path('api/job_convert_rds_to_h5ad_query/', JobConvertRdsToH5adQuery, name='api_job_convert_rds_to_h5ad_query'),
    re_path('api/job_convert_rds_to_h5ad_file_output_query/', JobConvertRdsToH5adFileOutputQuery, name='api_job_convert_rds_to_h5ad_file_output_query'),
    re_path('api/job_convert_rds_to_h5ad_file_output_query_by_id/', JobConvertRdsToH5adFileOutputQueryByID, name='api_job_convert_rds_to_h5ad_file_output_query_by_id'),
]
