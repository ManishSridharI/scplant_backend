import os
import datetime

from django.db import models
from django.core.validators import FileExtensionValidator
from django.dispatch import receiver

from accounts.models.CustomUserModel import CustomUserModel

from ..storage import OverwriteStorage


def job_annotate_and_plot_top25_markers_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(instance.job_annotate_and_plot_file_creation_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), filename)

def job_annotate_and_plot_marker_genes_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(instance.job_annotate_and_plot_file_creation_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), filename)

def job_annotate_and_plot_output_with_celltype_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(instance.job_annotate_and_plot_file_creation_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), filename)

def job_annotate_and_plot_prediction_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(instance.job_annotate_and_plot_file_creation_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), filename)

def job_annotate_and_plot_stdout_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/stdout_files/{2}'.format(instance.job_annotate_and_plot_file_creation_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), filename)

def job_annotate_and_plot_stderr_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/stderr_files/{2}'.format(instance.job_annotate_and_plot_file_creation_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")), filename)


class JobAnnotateAndPlotFileOutputModel(models.Model):
    job_annotate_and_plot_top25_markers_file = models.FileField(upload_to=job_annotate_and_plot_top25_markers_file_upload_directory, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_annotate_and_plot_marker_genes_file = models.FileField(upload_to=job_annotate_and_plot_marker_genes_file_upload_directory, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['csv'])])
    job_annotate_and_plot_output_with_celltype_file = models.FileField(upload_to=job_annotate_and_plot_output_with_celltype_file_upload_directory, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['h5ad'])])
    job_annotate_and_plot_prediction_file = models.FileField(upload_to=job_annotate_and_plot_prediction_file_upload_directory, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['csv'])])
    job_annotate_and_plot_stdout_file = models.FileField(upload_to=job_annotate_and_plot_stdout_file_upload_directory, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_annotate_and_plot_stderr_file = models.FileField(upload_to=job_annotate_and_plot_stderr_file_upload_directory, storage=OverwriteStorage(), validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_creation_timestamp = models.DateTimeField(auto_now_add=True)  # Automatically sets when created
    job_update_timestamp = models.DateTimeField(auto_now=True)  # Automatically updates on save
    job_annotate_and_plot_file_creation_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)


@receiver(models.signals.post_delete, sender=JobAnnotateAndPlotFileOutputModel)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    job_annotate_and_plot_top25_markers_file_path = instance.job_annotate_and_plot_top25_markers_file.path
    job_annotate_and_plot_marker_genes_file_path = instance.job_annotate_and_plot_marker_genes_file.path
    job_annotate_and_plot_output_with_celltype_file_path = instance.job_annotate_and_plot_output_with_celltype_file.path
    job_annotate_and_plot_prediction_file_path = instance.job_annotate_and_plot_prediction_file.path
    job_annotate_and_plot_stdout_file_path = instance.job_annotate_and_plot_stdout_file.path
    job_annotate_and_plot_stderr_file_path = instance.job_annotate_and_plot_stderr_file.path
    if instance.job_annotate_and_plot_top25_markers_file:
        instance.job_annotate_and_plot_top25_markers_file.delete(save=False)
    if instance.job_annotate_and_plot_marker_genes_file:
        instance.job_annotate_and_plot_marker_genes_file.delete(save=False)
    if instance.job_annotate_and_plot_output_with_celltype_file:
        instance.job_annotate_and_plot_output_with_celltype_file.delete(save=False)
    if instance.job_annotate_and_plot_prediction_file:
        instance.job_annotate_and_plot_prediction_file.delete(save=False)
    if instance.job_annotate_and_plot_stdout_file:
        instance.job_annotate_and_plot_stdout_file.delete(save=False)
    if instance.job_annotate_and_plot_stderr_file:
        instance.job_annotate_and_plot_stderr_file.delete(save=False)
    job_annotate_and_plot_folder_list = [
        os.path.dirname(os.path.dirname(job_annotate_and_plot_top25_markers_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_marker_genes_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_output_with_celltype_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_prediction_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_stdout_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_stderr_file_path))
    ]
    for i in range(0, len(job_annotate_and_plot_folder_list)):
        try:
            if os.path.exists(job_annotate_and_plot_folder_list[i]):
                os.rmdir(job_annotate_and_plot_folder_list[i])
        except Exception as e:
            pass
