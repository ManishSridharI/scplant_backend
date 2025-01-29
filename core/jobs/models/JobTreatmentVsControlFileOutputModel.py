import os
import datetime

from functools import partial

from django.db import models
from django.core.validators import FileExtensionValidator
from django.dispatch import receiver

from accounts.models.CustomUserModel import CustomUserModel
from organisms.models.OrganismModel import OrganismModel

from ..storage import OverwriteStorage


def job_treatment_vs_control_comp_celltype_dist_file_upload_directory(instance, filename):
    return 'jobs/treatment_vs_control/{0}/{1}/output/{2}'.format(
        instance.job_treatment_vs_control_file_creation_user.username,
        instance.job_treatment_vs_control_folder,
        filename
    )


def job_treatment_vs_control_control_vs_conditions_markers_file_upload_directory(instance, filename):
    return 'jobs/treatment_vs_control/{0}/{1}/output/{2}'.format(
        instance.job_treatment_vs_control_file_creation_user.username,
        instance.job_treatment_vs_control_folder,
        filename
    )


def job_treatment_vs_control_conditions_vs_control_markers_file_upload_directory(instance, filename):
    return 'jobs/treatment_vs_control/{0}/{1}/output/{2}'.format(
        instance.job_treatment_vs_control_file_creation_user.username,
        instance.job_treatment_vs_control_folder,
        filename
    )


def job_treatment_vs_control_stdout_file_upload_directory(instance, filename):
    return 'jobs/treatment_vs_control/{0}/{1}/stdout_files/{2}'.format(
        instance.job_treatment_vs_control_file_creation_user.username,
        instance.job_treatment_vs_control_folder,
        filename
    )


def job_treatment_vs_control_stderr_file_upload_directory(instance, filename):
    return 'jobs/treatment_vs_control/{0}/{1}/stderr_files/{2}'.format(
        instance.job_treatment_vs_control_file_creation_user.username,
        instance.job_treatment_vs_control_folder,
        filename
    )


def job_treatment_vs_control_control_vs_condition1_file_upload_directory(instance, filename):
    return 'jobs/treatment_vs_control/{0}/{1}/output/control_vs_condition1/{2}'.format(
        instance.job_treatment_vs_control_file_creation_user.username,
        instance.job_treatment_vs_control_folder,
        filename
    )


def job_treatment_vs_control_control_vs_condition2_file_upload_directory(instance, filename):
    return 'jobs/treatment_vs_control/{0}/{1}/output/control_vs_condition2/{2}'.format(
        instance.job_treatment_vs_control_file_creation_user.username,
        instance.job_treatment_vs_control_folder,
        filename
    )


class JobTreatmentVsControlFileOutputModel(models.Model):
    job_organism = models.ForeignKey(OrganismModel, on_delete=models.CASCADE, null=False)
    job_celery_task_id = models.CharField(max_length=255, unique=True, null=True, blank=False)
    job_treatment_vs_control_folder = models.CharField(max_length=200, null=True, blank=False)
    job_treatment_vs_control_folder_path = models.CharField(max_length=200, null=True, blank=False)

    job_treatment_vs_control_comp_celltype_dist_file = models.FileField(
        upload_to=job_treatment_vs_control_comp_celltype_dist_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['pdf'])],
        null=True
    )
    job_treatment_vs_control_control_vs_conditions_markers_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_conditions_markers_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )
    job_treatment_vs_control_conditions_vs_control_markers_file = models.FileField(
        upload_to=job_treatment_vs_control_conditions_vs_control_markers_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )

    job_treatment_vs_control_stdout_file = models.FileField(
        upload_to=job_treatment_vs_control_stdout_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['txt'])]
    )
    job_treatment_vs_control_stderr_file = models.FileField(
        upload_to=job_treatment_vs_control_stderr_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['txt'])]
    )

    job_treatment_vs_control_ctrl_vs_cond1_cond1_all_DEGs_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_condition1_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )
    job_treatment_vs_control_ctrl_vs_cond1_cond1_top10_DEGs_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_condition1_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )
    job_treatment_vs_control_ctrl_vs_cond1_cond1_top25_DEGs_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_condition1_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )
    job_treatment_vs_control_ctrl_vs_cond1_cond1_top5_DEGs_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_condition1_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )
    job_treatment_vs_control_ctrl_vs_cond1_ctrl_all_DEGs_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_condition1_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )
    job_treatment_vs_control_ctrl_vs_cond1_ctrl_top10_DEGs_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_condition1_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )
    job_treatment_vs_control_ctrl_vs_cond1_ctrl_top25_DEGs_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_condition1_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )
    job_treatment_vs_control_ctrl_vs_cond1_ctrl_top5_DEGs_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_condition1_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )

    job_treatment_vs_control_ctrl_vs_cond2_cond2_all_DEGs_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_condition2_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )
    job_treatment_vs_control_ctrl_vs_cond2_cond2_top10_DEGs_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_condition2_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )
    job_treatment_vs_control_ctrl_vs_cond2_cond2_top25_DEGs_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_condition2_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )
    job_treatment_vs_control_ctrl_vs_cond2_cond2_top5_DEGs_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_condition2_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )
    job_treatment_vs_control_ctrl_vs_cond2_ctrl_all_DEGs_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_condition2_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )
    job_treatment_vs_control_ctrl_vs_cond2_ctrl_top10_DEGs_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_condition2_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )
    job_treatment_vs_control_ctrl_vs_cond2_ctrl_top25_DEGs_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_condition2_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )
    job_treatment_vs_control_ctrl_vs_cond2_ctrl_top5_DEGs_file = models.FileField(
        upload_to=job_treatment_vs_control_control_vs_condition2_file_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['xlsx'])],
        null=True
    )

    job_treatment_vs_control_file_creation_timestamp = models.DateTimeField(auto_now_add=True)
    job_treatment_vs_control_file_creation_user = models.ForeignKey(
        CustomUserModel,
        on_delete=models.CASCADE,
        null=False
    )


@receiver(models.signals.post_save, sender=JobTreatmentVsControlFileOutputModel)
def your_function_name(sender, instance, created, **kwargs):
    if created:
        if instance.job_treatment_vs_control_folder:
            instance.job_treatment_vs_control_folder_path = os.path.join(
                'jobs/treatment_vs_control',
                instance.job_treatment_vs_control_file_creation_user.username,
                instance.job_treatment_vs_control_folder,
                'output'
            )
            instance.save()
    else:
        try:
            if instance.job_treatment_vs_control_folder:
                if os.path.join('jobs/treatment_vs_control', instance.job_treatment_vs_control_file_creation_user.username, instance.job_treatment_vs_control_folder, 'output') != instance.job_treatment_vs_control_folder_path:
                    instance.job_treatment_vs_control_folder_path = os.path.join(
                        'jobs/treatment_vs_control',
                        instance.job_treatment_vs_control_file_creation_user.username,
                        instance.job_treatment_vs_control_folder,
                        'output'
                    )
                    instance.save()
        except Exception as e:
            print(e)


@receiver(models.signals.post_delete, sender=JobTreatmentVsControlFileOutputModel)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    job_treatment_vs_control_comp_celltype_dist_file_path = instance.job_treatment_vs_control_comp_celltype_dist_file.path
    job_treatment_vs_control_control_vs_conditions_markers_file_path = instance.job_treatment_vs_control_control_vs_conditions_markers_file.path
    job_treatment_vs_control_conditions_vs_control_markers_file_path = instance.job_treatment_vs_control_conditions_vs_control_markers_file.path
    job_treatment_vs_control_stdout_file_path = instance.job_treatment_vs_control_stdout_file.path
    job_treatment_vs_control_stderr_file_path = instance.job_treatment_vs_control_stderr_file.path

    if instance.job_treatment_vs_control_comp_celltype_dist_file:
        instance.job_treatment_vs_control_comp_celltype_dist_file.delete(save=False)
    if instance.job_treatment_vs_control_control_vs_conditions_markers_file:
        instance.job_treatment_vs_control_control_vs_conditions_markers_file.delete(save=False)
    if instance.job_treatment_vs_control_conditions_vs_control_markers_file:
        instance.job_treatment_vs_control_conditions_vs_control_markers_file.delete(save=False)
    if instance.job_treatment_vs_control_stdout_file:
        instance.job_treatment_vs_control_stdout_file.delete(save=False)
    if instance.job_treatment_vs_control_stderr_file:
        instance.job_treatment_vs_control_stderr_file.delete(save=False)
    job_treatment_vs_control_folder_list = [
        os.path.dirname(os.path.dirname(job_treatment_vs_control_comp_celltype_dist_file_path)),
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
