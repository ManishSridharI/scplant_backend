import json

from django.db import models

from accounts.models.CustomUserModel import CustomUserModel
from scripts.models.ScriptModel import ScriptModel
from datasets.models.DatasetModel import DatasetModel
from predictors.models.PredictorModel import PredictorModel

from .JobAnnotateAndPlotArgumentModel import JobAnnotateAndPlotArgumentModel
from .JobAnnotateAndPlotFileOutputModel import JobAnnotateAndPlotFileOutputModel


class JobAnnotateAndPlotModel(models.Model):
    job_name = models.CharField(max_length=200, null=False, blank=False)
    job_script = models.ForeignKey(ScriptModel, on_delete=models.CASCADE, null=False)
    job_dataset = models.ForeignKey(DatasetModel, on_delete=models.CASCADE, null=False)
    job_predictor = models.ForeignKey(PredictorModel, on_delete=models.CASCADE, null=False)
    job_annotate_and_plot_argument = models.ForeignKey(JobAnnotateAndPlotArgumentModel, on_delete=models.CASCADE, null=False)
    job_annotate_and_plot_file_output = models.ForeignKey(JobAnnotateAndPlotFileOutputModel, on_delete=models.CASCADE, null=False)
    job_celery_task_id = models.CharField(max_length=255, unique=True, null=False, blank=False)
    job_celery_task_status = models.CharField(max_length=50, null=True, blank=False, db_index=True)
    job_celery_task_result = models.TextField(null=True, blank=False)
    job_creation_timestamp = models.DateTimeField(auto_now_add=True)  # Automatically sets when created
    job_creation_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)

    def __str__(self):
        output_string = json.dumps({
            "job_name": self.job_name,
            "job_celery_task_id": self.job_celery_task_id,
            "job_celery_task_status": self.job_celery_task_status,
            "job_celery_task_result": self.job_celery_task_result
        })
        return output_string

    @property
    def get_job_name(self):
        return f"{self.job_name}"
