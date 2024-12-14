import json

from django.db import models

from accounts.models.CustomUserModel import CustomUserModel
from datasets.models.DatasetModel import DatasetModel
from predictors.models.PredictorModel import PredictorModel


class JobModel(models.Model):
    job_name = models.CharField(max_length=200, null=False, blank=False)
    job_dataset = models.ForeignKey(DatasetModel, on_delete=models.CASCADE, null=True)
    job_predictor = models.ForeignKey(PredictorModel, on_delete=models.CASCADE, null=True)
    job_celery_task = models.CharField(max_length=255, unique=True, null=False, blank=False)
    job_creation_timestamp = models.DateTimeField(auto_now_add=True)  # Automatically sets when created
    job_creation_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)

    def __str__(self):
        output_string = json.dumps({
            "job_name": self.job_name,
        })
        return output_string

    @property
    def get_job_name(self):
        return f"{self.job_name}"
