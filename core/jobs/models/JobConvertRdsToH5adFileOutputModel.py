import os
import datetime

from functools import partial

from django.db import models
from django.core.validators import FileExtensionValidator
from django.dispatch import receiver

from accounts.models.CustomUserModel import CustomUserModel
from organisms.models.OrganismModel import OrganismModel

from ..storage import OverwriteStorage


def job_convert_rds_to_h5ad_output_file_upload_directory(instance, filename):
    return 'jobs/convert_rds_to_h5ad/{0}/{1}/output_files/{2}'.format(instance.job_convert_rds_to_h5ad_file_creation_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), filename)

def job_convert_rds_to_h5ad_stdout_file_upload_directory(instance, filename):
    return 'jobs/convert_rds_to_h5ad/{0}/{1}/stdout_files/{2}'.format(instance.job_convert_rds_to_h5ad_file_creation_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), filename)

def job_convert_rds_to_h5ad_stderr_file_upload_directory(instance, filename):
    return 'jobs/convert_rds_to_h5ad/{0}/{1}/stderr_files/{2}'.format(instance.job_convert_rds_to_h5ad_file_creation_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), filename)


class JobConvertRdsToH5adFileOutputModel(models.Model):
    job_organism = models.ForeignKey(OrganismModel, on_delete=models.CASCADE, null=False)
    job_celery_task_id = models.CharField(max_length=255, unique=True, null=True, blank=False)
    job_convert_rds_to_h5ad_output_file = models.FileField(upload_to=job_convert_rds_to_h5ad_output_file_upload_directory, max_length=200, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['h5ad'])])
    job_convert_rds_to_h5ad_stdout_file = models.FileField(upload_to=job_convert_rds_to_h5ad_stdout_file_upload_directory, max_length=200, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_convert_rds_to_h5ad_stderr_file = models.FileField(upload_to=job_convert_rds_to_h5ad_stderr_file_upload_directory, max_length=200, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_convert_rds_to_h5ad_file_creation_timestamp = models.DateTimeField(auto_now_add=True)
    job_convert_rds_to_h5ad_file_creation_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)


@receiver(models.signals.post_delete, sender=JobConvertRdsToH5adFileOutputModel)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    job_convert_rds_to_h5ad_output_file_path = instance.job_convert_rds_to_h5ad_output_file.path
    job_convert_rds_to_h5ad_stdout_file_path = instance.job_convert_rds_to_h5ad_stdout_file.path
    job_convert_rds_to_h5ad_stderr_file_path = instance.job_convert_rds_to_h5ad_stderr_file.path
    if instance.job_convert_rds_to_h5ad_output_file:
        instance.job_convert_rds_to_h5ad_output_file.delete(save=False)
    if instance.job_convert_rds_to_h5ad_stdout_file:
        instance.job_convert_rds_to_h5ad_stdout_file.delete(save=False)
    if instance.job_convert_rds_to_h5ad_stderr_file:
        instance.job_convert_rds_to_h5ad_stderr_file.delete(save=False)
    job_convert_rds_to_h5ad_folder_list = [
        os.path.dirname(os.path.dirname(job_convert_rds_to_h5ad_output_file_path)),
        os.path.dirname(os.path.dirname(job_convert_rds_to_h5ad_stdout_file_path)),
        os.path.dirname(os.path.dirname(job_convert_rds_to_h5ad_stderr_file_path))
    ]
    for i in range(0, len(job_convert_rds_to_h5ad_folder_list)):
        try:
            if os.path.exists(job_convert_rds_to_h5ad_folder_list[i]):
                os.rmdir(job_convert_rds_to_h5ad_folder_list[i])
        except Exception as e:
            pass
