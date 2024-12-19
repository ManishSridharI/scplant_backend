import os
import datetime

from functools import partial

from django.db import models
from django.core.validators import FileExtensionValidator
from django.dispatch import receiver

from accounts.models.CustomUserModel import CustomUserModel

from ..storage import OverwriteStorage


def job_treatment_vs_control_condition1_marker_genes_file_upload_directory(instance, filename, now_strftime):
    return 'jobs/treatment_vs_control/{0}/{1}/output/{2}'.format(instance.job_treatment_vs_control_file_creation_user.username, now_strftime, filename)

def job_treatment_vs_control_condition1_top25_markers_file_upload_directory(instance, filename, now_strftime):
    return 'jobs/treatment_vs_control/{0}/{1}/output/{2}'.format(instance.job_treatment_vs_control_file_creation_user.username, now_strftime, filename)

def job_treatment_vs_control_condition1_top10_genes_dotplot_file_upload_directory(instance, filename, now_strftime):
    return 'jobs/treatment_vs_control/{0}/{1}/output/{2}'.format(instance.job_treatment_vs_control_file_creation_user.username, now_strftime, filename)

def job_treatment_vs_control_condition2_marker_genes_file_upload_directory(instance, filename, now_strftime):
    return 'jobs/treatment_vs_control/{0}/{1}/output/{2}'.format(instance.job_treatment_vs_control_file_creation_user.username, now_strftime, filename)

def job_treatment_vs_control_condition2_top25_markers_file_upload_directory(instance, filename, now_strftime):
    return 'jobs/treatment_vs_control/{0}/{1}/output/{2}'.format(instance.job_treatment_vs_control_file_creation_user.username, now_strftime, filename)

def job_treatment_vs_control_condition2_top10_genes_dotplot_file_upload_directory(instance, filename, now_strftime):
    return 'jobs/treatment_vs_control/{0}/{1}/output/{2}'.format(instance.job_treatment_vs_control_file_creation_user.username, now_strftime, filename)

def job_treatment_vs_control_control_vs_conditions_markers_file_upload_directory(instance, filename, now_strftime):
    return 'jobs/treatment_vs_control/{0}/{1}/output/{2}'.format(instance.job_treatment_vs_control_file_creation_user.username, now_strftime, filename)

def job_treatment_vs_control_conditions_vs_control_markers_file_upload_directory(instance, filename, now_strftime):
    return 'jobs/treatment_vs_control/{0}/{1}/output/{2}'.format(instance.job_treatment_vs_control_file_creation_user.username, now_strftime, filename)

def job_treatment_vs_control_stdout_file_upload_directory(instance, filename, now_strftime):
    return 'jobs/treatment_vs_control/{0}/{1}/stdout_files/{2}'.format(instance.job_treatment_vs_control_file_creation_user.username, now_strftime, filename)

def job_treatment_vs_control_stderr_file_upload_directory(instance, filename, now_strftime):
    return 'jobs/treatment_vs_control/{0}/{1}/stderr_files/{2}'.format(instance.job_treatment_vs_control_file_creation_user.username, now_strftime, filename)


class JobTreatmentVsControlFileOutputModel(models.Model):
    _now_strftime = str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S"))
    job_celery_task_id = models.CharField(max_length=255, unique=True, null=True, blank=False)
    job_treatment_vs_control_condition1_marker_genes_file = models.FileField(upload_to=partial(job_treatment_vs_control_condition1_marker_genes_file_upload_directory, now_strftime=_now_strftime), max_length=200, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['csv'])])
    job_treatment_vs_control_condition1_top25_markers_file = models.FileField(upload_to=partial(job_treatment_vs_control_condition1_top25_markers_file_upload_directory, now_strftime=_now_strftime), max_length=200, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_treatment_vs_control_condition1_top10_genes_dotplot_file = models.FileField(upload_to=partial(job_treatment_vs_control_condition1_top10_genes_dotplot_file_upload_directory, now_strftime=_now_strftime), max_length=200, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['pdf'])])
    job_treatment_vs_control_condition2_marker_genes_file = models.FileField(upload_to=partial(job_treatment_vs_control_condition2_marker_genes_file_upload_directory, now_strftime=_now_strftime), max_length=200, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['csv'])])
    job_treatment_vs_control_condition2_top25_markers_file = models.FileField(upload_to=partial(job_treatment_vs_control_condition2_top25_markers_file_upload_directory, now_strftime=_now_strftime), max_length=200, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_treatment_vs_control_condition2_top10_genes_dotplot_file = models.FileField(upload_to=partial(job_treatment_vs_control_condition2_top10_genes_dotplot_file_upload_directory, now_strftime=_now_strftime), max_length=200, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['pdf'])])
    job_treatment_vs_control_control_vs_conditions_markers_file = models.FileField(upload_to=partial(job_treatment_vs_control_control_vs_conditions_markers_file_upload_directory, now_strftime=_now_strftime), max_length=200, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_treatment_vs_control_conditions_vs_control_markers_file = models.FileField(upload_to=partial(job_treatment_vs_control_conditions_vs_control_markers_file_upload_directory, now_strftime=_now_strftime), max_length=200, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_treatment_vs_control_stdout_file = models.FileField(upload_to=partial(job_treatment_vs_control_stdout_file_upload_directory, now_strftime=_now_strftime), max_length=200, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_treatment_vs_control_stderr_file = models.FileField(upload_to=partial(job_treatment_vs_control_stderr_file_upload_directory, now_strftime=_now_strftime), max_length=200, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_treatment_vs_control_file_creation_timestamp = models.DateTimeField(auto_now_add=True)
    job_treatment_vs_control_file_creation_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)


@receiver(models.signals.post_delete, sender=JobTreatmentVsControlFileOutputModel)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    job_treatment_vs_control_condition1_marker_genes_file_path = instance.job_treatment_vs_control_condition1_marker_genes_file.path
    job_treatment_vs_control_condition1_top25_markers_file_path = instance.job_treatment_vs_control_condition1_top25_markers_file.path
    job_treatment_vs_control_condition1_top10_genes_dotplot_file_path = instance.job_treatment_vs_control_condition1_top10_genes_dotplot_file.path
    job_treatment_vs_control_condition2_marker_genes_file_path = instance.job_treatment_vs_control_condition2_marker_genes_file.path
    job_treatment_vs_control_condition2_top25_markers_file_path = instance.job_treatment_vs_control_condition2_top25_markers_file.path
    job_treatment_vs_control_condition2_top10_genes_dotplot_file_path = instance.job_treatment_vs_control_condition2_top10_genes_dotplot_file.path
    job_treatment_vs_control_control_vs_conditions_markers_file_path = instance.job_treatment_vs_control_control_vs_conditions_markers_file.path
    job_treatment_vs_control_conditions_vs_control_markers_file_path = instance.job_treatment_vs_control_conditions_vs_control_markers_file.path
    job_treatment_vs_control_stdout_file_path = instance.job_treatment_vs_control_stdout_file.path
    job_treatment_vs_control_stderr_file_path = instance.job_treatment_vs_control_stderr_file.path
    if instance.job_treatment_vs_control_condition1_marker_genes_file:
        instance.job_treatment_vs_control_condition1_marker_genes_file.delete(save=False)
    if instance.job_treatment_vs_control_condition1_top25_markers_file:
        instance.job_treatment_vs_control_condition1_top25_markers_file.delete(save=False)
    if instance.job_treatment_vs_control_condition1_top10_genes_dotplot_file:
        instance.job_treatment_vs_control_condition1_top10_genes_dotplot_file.delete(save=False)
    if instance.job_treatment_vs_control_condition2_marker_genes_file:
        instance.job_treatment_vs_control_condition2_marker_genes_file.delete(save=False)
    if instance.job_treatment_vs_control_condition2_top25_markers_file:
        instance.job_treatment_vs_control_condition2_top25_markers_file.delete(save=False)
    if instance.job_treatment_vs_control_condition2_top10_genes_dotplot_file:
        instance.job_treatment_vs_control_condition2_top10_genes_dotplot_file.delete(save=False)
    if instance.job_treatment_vs_control_control_vs_conditions_markers_file:
        instance.job_treatment_vs_control_control_vs_conditions_markers_file.delete(save=False)
    if instance.job_treatment_vs_control_conditions_vs_control_markers_file:
        instance.job_treatment_vs_control_conditions_vs_control_markers_file.delete(save=False)
    if instance.job_treatment_vs_control_stdout_file:
        instance.job_treatment_vs_control_stdout_file.delete(save=False)
    if instance.job_treatment_vs_control_stderr_file:
        instance.job_treatment_vs_control_stderr_file.delete(save=False)
    job_treatment_vs_control_folder_list = [
        os.path.dirname(os.path.dirname(job_treatment_vs_control_condition1_marker_genes_file_path)),
        os.path.dirname(os.path.dirname(job_treatment_vs_control_condition1_top25_markers_file_path)),
        os.path.dirname(os.path.dirname(job_treatment_vs_control_condition1_top10_genes_dotplot_file_path)),
        os.path.dirname(os.path.dirname(job_treatment_vs_control_condition2_marker_genes_file_path)),
        os.path.dirname(os.path.dirname(job_treatment_vs_control_condition2_top25_markers_file_path)),
        os.path.dirname(os.path.dirname(job_treatment_vs_control_condition2_top10_genes_dotplot_file_path)),
        os.path.dirname(os.path.dirname(job_treatment_vs_control_control_vs_conditions_markers_file_path)),
        os.path.dirname(os.path.dirname(job_treatment_vs_control_conditions_vs_control_markers_file_path)),
        os.path.dirname(os.path.dirname(job_treatment_vs_control_stdout_file_path)),
        os.path.dirname(os.path.dirname(job_treatment_vs_control_stderr_file_path))
    ]
    for i in range(0, len(job_treatment_vs_control_folder_list)):
        try:
            if os.path.exists(job_treatment_vs_control_folder_list[i]):
                os.rmdir(job_treatment_vs_control_folder_list[i])
        except Exception as e:
            pass
