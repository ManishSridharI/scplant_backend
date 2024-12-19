import os
import datetime

from functools import partial

from django.db import models
from django.core.validators import FileExtensionValidator
from django.dispatch import receiver

from accounts.models.CustomUserModel import CustomUserModel

from ..storage import OverwriteStorage


def job_compare_cell_type_dist_output_file_upload_directory(instance, filename):
    return 'jobs/compare_cell_type_dist/{0}/{1}/output/{2}'.format(instance.job_compare_cell_type_dist_file_creation_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), filename)

def job_compare_cell_type_dist_stdout_file_upload_directory(instance, filename):
    return 'jobs/compare_cell_type_dist/{0}/{1}/stdout_files/{2}'.format(instance.job_compare_cell_type_dist_file_creation_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), filename)

def job_compare_cell_type_dist_stderr_file_upload_directory(instance, filename):
    return 'jobs/compare_cell_type_dist/{0}/{1}/stderr_files/{2}'.format(instance.job_compare_cell_type_dist_file_creation_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), filename)


class JobCompareCellTypeDistFileOutputModel(models.Model):
    job_celery_task_id = models.CharField(max_length=255, unique=True, null=True, blank=False)
    job_compare_cell_type_dist_output_file = models.FileField(upload_to=job_compare_cell_type_dist_output_file_upload_directory,max_length=200, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['pdf'])], null=True)
    job_compare_cell_type_dist_stdout_file = models.FileField(upload_to=job_compare_cell_type_dist_stdout_file_upload_directory, max_length=200, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_compare_cell_type_dist_stderr_file = models.FileField(upload_to=job_compare_cell_type_dist_stderr_file_upload_directory, max_length=200, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_compare_cell_type_dist_file_creation_timestamp = models.DateTimeField(auto_now_add=True)
    job_compare_cell_type_dist_file_creation_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)


@receiver(models.signals.post_delete, sender=JobCompareCellTypeDistFileOutputModel)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    job_compare_cell_type_dist_output_file_path = instance.job_compare_cell_type_dist_output_file.path
    job_compare_cell_type_dist_stdout_file_path = instance.job_compare_cell_type_dist_stdout_file.path
    job_compare_cell_type_dist_stderr_file_path = instance.job_compare_cell_type_dist_stderr_file.path
    if instance.job_compare_cell_type_dist_output_file:
        instance.job_compare_cell_type_dist_output_file.delete(save=False)
    if instance.job_compare_cell_type_dist_stdout_file:
        instance.job_compare_cell_type_dist_stdout_file.delete(save=False)
    if instance.job_compare_cell_type_dist_stderr_file:
        instance.job_compare_cell_type_dist_stderr_file.delete(save=False)
    job_compare_cell_type_dist_folder_list = [
        os.path.dirname(os.path.dirname(
            job_compare_cell_type_dist_output_file_path)),
        os.path.dirname(os.path.dirname(
            job_compare_cell_type_dist_stdout_file_path)),
        os.path.dirname(os.path.dirname(
            job_compare_cell_type_dist_stderr_file_path))
    ]
    for i in range(0, len(job_compare_cell_type_dist_folder_list)):
        try:
            if os.path.exists(job_compare_cell_type_dist_folder_list[i]):
                os.rmdir(job_compare_cell_type_dist_folder_list[i])
        except Exception as e:
            pass
