from django.db import models

from accounts.models.CustomUserModel import CustomUserModel


class JobAnnotateAndPlotArgumentModel(models.Model):
    job_celery_task_id = models.CharField(max_length=255, unique=True, null=True, blank=False)
    job_annotate_and_plot_gene_number = models.IntegerField(null=False, blank=False)
    job_annotate_and_plot_argument_creation_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)
