from django.db import models
from django.core.validators import FileExtensionValidator
from django.dispatch import receiver

from accounts.models.CustomUserModel import CustomUserModel


def job_inference_log_file_upload_directory(instance, filename):
    return 'jobs/inference/log_files/{0}/{1}'.format(instance.job_inference_file_creation_user.username, filename)

def job_inference_prediction_file_upload_directory(instance, filename):
    return 'jobs/inference/prediction_files/{0}/{1}'.format(instance.job_inference_file_creation_user.username, filename)

def job_inference_stdout_file_upload_directory(instance, filename):
    return 'jobs/inference/stdout_files/{0}/{1}'.format(instance.job_inference_file_creation_user.username, filename)

def job_inference_stderr_file_upload_directory(instance, filename):
    return 'jobs/inference/stderr_files/{0}/{1}'.format(instance.job_inference_file_creation_user.username, filename)


class JobInferenceFileOutputModel(models.Model):
    job_inference_log_file = models.FileField(upload_to=job_inference_log_file_upload_directory, validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_inference_prediction_file = models.FileField(upload_to=job_inference_prediction_file_upload_directory, validators=[FileExtensionValidator(allowed_extensions=['csv'])])
    job_inference_stdout_file = models.FileField(upload_to=job_inference_stdout_file_upload_directory, validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_inference_stderr_file = models.FileField(upload_to=job_inference_stderr_file_upload_directory, validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_creation_timestamp = models.DateTimeField(auto_now_add=True)  # Automatically sets when created
    job_update_timestamp = models.DateTimeField(auto_now=True)  # Automatically updates on save
    job_inference_file_creation_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)



@receiver(models.signals.post_delete, sender=JobInferenceFileOutputModel)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    if instance.job_inference_log_file:
        instance.job_inference_log_file.delete(save=False)
    if instance.job_inference_prediction_file:
        instance.job_inference_prediction_file.delete(save=False)
    if instance.job_inference_stdout_file:
        instance.job_inference_stdout_file.delete(save=False)
    if instance.job_inference_stderr_file:
        instance.job_inference_stderr_file.delete(save=False)
