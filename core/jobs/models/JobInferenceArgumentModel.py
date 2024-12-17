from django.db import models

from accounts.models.CustomUserModel import CustomUserModel


class JobInferenceArgumentModel(models.Model):
    job_inference_gene_number = models.IntegerField(null=False, blank=False)
    job_inference_argument_creation_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)
