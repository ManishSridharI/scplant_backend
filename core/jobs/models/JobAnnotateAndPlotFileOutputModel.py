import os
import shutil
import datetime

from functools import partial

from django.db import models
from django.core.validators import FileExtensionValidator
from django.dispatch import receiver

from accounts.models.CustomUserModel import CustomUserModel
from organisms.models.OrganismModel import OrganismModel

from ..storage import OverwriteStorage


def job_annotate_and_plot_output_with_celltype_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(
        instance.job_annotate_and_plot_file_creation_user.username,
        instance.job_annotate_and_plot_folder,
        filename
    )


def job_annotate_and_plot_prediction_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(
        instance.job_annotate_and_plot_file_creation_user.username,
        instance.job_annotate_and_plot_folder,
        filename
    )


def job_annotate_and_plot_log_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(
        instance.job_annotate_and_plot_file_creation_user.username,
        instance.job_annotate_and_plot_folder,
        filename
    )


def job_annotate_and_plot_stats_csv_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(
        instance.job_annotate_and_plot_file_creation_user.username,
        instance.job_annotate_and_plot_folder,
        filename
    )


def job_annotate_and_plot_stats_pdf_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(
        instance.job_annotate_and_plot_file_creation_user.username,
        instance.job_annotate_and_plot_folder,
        filename
    )


def job_annotate_and_plot_marker_genes_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(
        instance.job_annotate_and_plot_file_creation_user.username,
        instance.job_annotate_and_plot_folder,
        filename
    )


def job_annotate_and_plot_annotate_tsne_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(
        instance.job_annotate_and_plot_file_creation_user.username,
        instance.job_annotate_and_plot_folder,
        filename
    )


def job_annotate_and_plot_annotate_umap_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(
        instance.job_annotate_and_plot_file_creation_user.username,
        instance.job_annotate_and_plot_folder,
        filename
    )


def job_annotate_and_plot_top3_genes_dotplot_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(
        instance.job_annotate_and_plot_file_creation_user.username,
        instance.job_annotate_and_plot_folder,
        filename
    )


def job_annotate_and_plot_all_markers_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(
        instance.job_annotate_and_plot_file_creation_user.username,
        instance.job_annotate_and_plot_folder,
        filename
    )


def job_annotate_and_plot_top5_markers_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(
        instance.job_annotate_and_plot_file_creation_user.username,
        instance.job_annotate_and_plot_folder,
        filename
    )


def job_annotate_and_plot_top10_markers_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(
        instance.job_annotate_and_plot_file_creation_user.username,
        instance.job_annotate_and_plot_folder,
        filename
    )


def job_annotate_and_plot_top25_markers_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/output/{2}'.format(
        instance.job_annotate_and_plot_file_creation_user.username,
        instance.job_annotate_and_plot_folder,
        filename
    )


def job_annotate_and_plot_stdout_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/stdout_files/{2}'.format(
        instance.job_annotate_and_plot_file_creation_user.username,
        instance.job_annotate_and_plot_folder,
        filename
    )


def job_annotate_and_plot_stderr_file_upload_directory(instance, filename):
    return 'jobs/annotate_and_plot/{0}/{1}/stderr_files/{2}'.format(
        instance.job_annotate_and_plot_file_creation_user.username,
        instance.job_annotate_and_plot_folder,
        filename
    )


class JobAnnotateAndPlotFileOutputModel(models.Model):
    job_organism = models.ForeignKey(OrganismModel, on_delete=models.CASCADE, null=False)
    job_celery_task_id = models.CharField(max_length=255, unique=True, null=True, blank=False)
    job_annotate_and_plot_folder = models.CharField(max_length=200, null=True, blank=False)
    job_annotate_and_plot_folder_path = models.CharField(max_length=200, null=True, blank=False)
    job_annotate_and_plot_output_with_celltype_file = models.FileField(
        upload_to=job_annotate_and_plot_output_with_celltype_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['h5ad'])])
    job_annotate_and_plot_prediction_file = models.FileField(
        upload_to=job_annotate_and_plot_prediction_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['csv'])])
    job_annotate_and_plot_log_file = models.FileField(
        upload_to=job_annotate_and_plot_log_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_annotate_and_plot_stats_csv_file = models.FileField(
        upload_to=job_annotate_and_plot_stats_csv_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['csv'])])
    job_annotate_and_plot_stats_pdf_file = models.FileField(
        upload_to=job_annotate_and_plot_stats_pdf_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['pdf'])])
    job_annotate_and_plot_marker_genes_file = models.FileField(
        upload_to=job_annotate_and_plot_marker_genes_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['csv'])])
    job_annotate_and_plot_annotate_tsne_file = models.FileField(
        upload_to=job_annotate_and_plot_annotate_tsne_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['pdf'])])
    job_annotate_and_plot_annotate_umap_file = models.FileField(
        upload_to=job_annotate_and_plot_annotate_umap_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['pdf'])])
    job_annotate_and_plot_top3_genes_dotplot_file = models.FileField(
        upload_to=job_annotate_and_plot_top3_genes_dotplot_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['pdf'])])
    job_annotate_and_plot_all_markers_file = models.FileField(
        upload_to=job_annotate_and_plot_all_markers_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])])
    job_annotate_and_plot_top5_markers_file = models.FileField(
        upload_to=job_annotate_and_plot_top5_markers_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])])
    job_annotate_and_plot_top10_markers_file = models.FileField(
        upload_to=job_annotate_and_plot_top10_markers_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])])
    job_annotate_and_plot_top25_markers_file = models.FileField(
        upload_to=job_annotate_and_plot_top25_markers_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])])
    job_annotate_and_plot_stdout_file = models.FileField(
        upload_to=job_annotate_and_plot_stdout_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_annotate_and_plot_stderr_file = models.FileField(
        upload_to=job_annotate_and_plot_stderr_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    job_annotate_and_plot_file_creation_timestamp = models.DateTimeField(auto_now_add=True)
    job_annotate_and_plot_file_creation_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)


@receiver(models.signals.post_save, sender=JobAnnotateAndPlotFileOutputModel)
def your_function_name(sender, instance, created, **kwargs):
    if created:
        if instance.job_annotate_and_plot_folder:
            instance.job_annotate_and_plot_folder_path = os.path.join(
                'jobs/annotate_and_plot',
                instance.job_annotate_and_plot_file_creation_user.username,
                instance.job_annotate_and_plot_folder,
                'output'
            )
            instance.save()
    else:
        try:
            if instance.job_annotate_and_plot_folder:
                if os.path.join('jobs/annotate_and_plot', instance.job_annotate_and_plot_file_creation_user.username, instance.job_annotate_and_plot_folder, 'output') != instance.job_annotate_and_plot_folder_path:
                    instance.job_annotate_and_plot_folder_path = os.path.join(
                        'jobs/annotate_and_plot',
                        instance.job_annotate_and_plot_file_creation_user.username,
                        instance.job_annotate_and_plot_folder,
                        'output'
                    )
                    instance.save()
        except Exception as e:
            print(e)


@receiver(models.signals.post_delete, sender=JobAnnotateAndPlotFileOutputModel)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    job_annotate_and_plot_output_with_celltype_file_path = instance.job_annotate_and_plot_output_with_celltype_file.path
    job_annotate_and_plot_prediction_file_path = instance.job_annotate_and_plot_prediction_file.path
    job_annotate_and_plot_log_file_path = instance.job_annotate_and_plot_log_file.path
    job_annotate_and_plot_stats_csv_file_path = instance.job_annotate_and_plot_stats_csv_file.path
    job_annotate_and_plot_stats_pdf_file_path = instance.job_annotate_and_plot_stats_pdf_file.path
    job_annotate_and_plot_marker_genes_file_path = instance.job_annotate_and_plot_marker_genes_file.path
    job_annotate_and_plot_annotate_tsne_file_path = instance.job_annotate_and_plot_annotate_tsne_file.path
    job_annotate_and_plot_annotate_umap_file_path = instance.job_annotate_and_plot_annotate_umap_file.path
    job_annotate_and_plot_top3_genes_dotplot_file_path = instance.job_annotate_and_plot_top3_genes_dotplot_file.path
    job_annotate_and_plot_all_markers_file_path = instance.job_annotate_and_plot_all_markers_file.path
    job_annotate_and_plot_top5_markers_file_path = instance.job_annotate_and_plot_top5_markers_file.path
    job_annotate_and_plot_top10_markers_file_path = instance.job_annotate_and_plot_top10_markers_file.path
    job_annotate_and_plot_top25_markers_file_path = instance.job_annotate_and_plot_top25_markers_file.path
    job_annotate_and_plot_stdout_file_path = instance.job_annotate_and_plot_stdout_file.path
    job_annotate_and_plot_stderr_file_path = instance.job_annotate_and_plot_stderr_file.path

    if instance.job_annotate_and_plot_output_with_celltype_file:
        instance.job_annotate_and_plot_output_with_celltype_file.delete(save=False)
    if instance.job_annotate_and_plot_prediction_file:
        instance.job_annotate_and_plot_prediction_file.delete(save=False)
    if instance.job_annotate_and_plot_log_file:
        instance.job_annotate_and_plot_log_file.delete(save=False)
    if instance.job_annotate_and_plot_stats_csv_file:
        instance.job_annotate_and_plot_stats_csv_file.delete(save=False)
    if instance.job_annotate_and_plot_stats_pdf_file:
        instance.job_annotate_and_plot_stats_pdf_file.delete(save=False)
    if instance.job_annotate_and_plot_marker_genes_file:
        instance.job_annotate_and_plot_marker_genes_file.delete(save=False)
    if instance.job_annotate_and_plot_annotate_tsne_file:
        instance.job_annotate_and_plot_annotate_tsne_file.delete(save=False)
    if instance.job_annotate_and_plot_annotate_umap_file:
        instance.job_annotate_and_plot_annotate_umap_file.delete(save=False)
    if instance.job_annotate_and_plot_top3_genes_dotplot_file:
        instance.job_annotate_and_plot_top3_genes_dotplot_file.delete(save=False)
    if instance.job_annotate_and_plot_all_markers_file:
        instance.job_annotate_and_plot_all_markers_file.delete(save=False)
    if instance.job_annotate_and_plot_top5_markers_file:
        instance.job_annotate_and_plot_top5_markers_file.delete(save=False)
    if instance.job_annotate_and_plot_top10_markers_file:
        instance.job_annotate_and_plot_top10_markers_file.delete(save=False)
    if instance.job_annotate_and_plot_top25_markers_file:
        instance.job_annotate_and_plot_top25_markers_file.delete(save=False)
    if instance.job_annotate_and_plot_stdout_file:
        instance.job_annotate_and_plot_stdout_file.delete(save=False)
    if instance.job_annotate_and_plot_stderr_file:
        instance.job_annotate_and_plot_stderr_file.delete(save=False)
    job_annotate_and_plot_folder_list = [
        os.path.dirname(os.path.dirname(job_annotate_and_plot_output_with_celltype_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_prediction_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_log_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_stats_csv_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_stats_pdf_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_marker_genes_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_annotate_tsne_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_annotate_umap_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_top3_genes_dotplot_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_all_markers_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_top5_markers_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_top10_markers_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_top25_markers_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_stdout_file_path)),
        os.path.dirname(os.path.dirname(job_annotate_and_plot_stderr_file_path))
    ]
    for i in range(0, len(job_annotate_and_plot_folder_list)):
        try:
            if os.path.exists(job_annotate_and_plot_folder_list[i]):
                shutil.rmtree(job_annotate_and_plot_folder_list[i])
        except Exception as e:
            pass
