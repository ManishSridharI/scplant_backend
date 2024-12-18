import os
import datetime

from django.db import models
from django.core.validators import FileExtensionValidator
from django.dispatch import receiver

from accounts.models.CustomUserModel import CustomUserModel

from ..storage import OverwriteStorage


def job_inference_log_file_upload_directory(instance, filename):
    return 'jobs/inference/{0}/{1}/log_files/{2}'.format(instance.job_inference_file_creation_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), filename)

def job_inference_prediction_file_upload_directory(instance, filename):
    return 'jobs/inference/{0}/{1}/prediction_files/{2}'.format(instance.job_inference_file_creation_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), filename)

def job_inference_stdout_file_upload_directory(instance, filename):
    return 'jobs/inference/{0}/{1}/stdout_files/{2}'.format(instance.job_inference_file_creation_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), filename)

def job_inference_stderr_file_upload_directory(instance, filename):
    return 'jobs/inference/{0}/{1}/stderr_files/{2}'.format(instance.job_inference_file_creation_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), filename)


class JobInferenceFileOutputModel(models.Model):
    job_celery_task_id = models.CharField(max_length=255, unique=True, null=True, blank=False)
    job_inference_log_file = models.FileField(upload_to=job_inference_log_file_upload_directory, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_inference_prediction_file = models.FileField(upload_to=job_inference_prediction_file_upload_directory, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['csv'])])
    job_inference_stdout_file = models.FileField(upload_to=job_inference_stdout_file_upload_directory, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_inference_stderr_file = models.FileField(upload_to=job_inference_stderr_file_upload_directory, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_inference_file_creation_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)


@receiver(models.signals.post_delete, sender=JobInferenceFileOutputModel)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    job_inference_log_file_path = instance.job_inference_log_file.path
    job_inference_prediction_file_path = instance.job_inference_prediction_file.path
    job_inference_stdout_file_path = instance.job_inference_stdout_file.path
    job_inference_stderr_file_path = instance.job_inference_stderr_file.path
    if instance.job_inference_log_file:
        instance.job_inference_log_file.delete(save=False)
    if instance.job_inference_prediction_file:
        instance.job_inference_prediction_file.delete(save=False)
    if instance.job_inference_stdout_file:
        instance.job_inference_stdout_file.delete(save=False)
    if instance.job_inference_stderr_file:
        instance.job_inference_stderr_file.delete(save=False)
    job_inference_folder_list = [
        os.path.dirname(os.path.dirname(job_inference_log_file_path)),
        os.path.dirname(os.path.dirname(job_inference_prediction_file_path)),
        os.path.dirname(os.path.dirname(job_inference_stdout_file_path)),
        os.path.dirname(os.path.dirname(job_inference_stderr_file_path))
    ]
    for i in range(0, len(job_inference_folder_list)):
        try:
            if os.path.exists(job_inference_folder_list[i]):
                os.rmdir(job_inference_folder_list[i])
        except Exception as e:
            pass
