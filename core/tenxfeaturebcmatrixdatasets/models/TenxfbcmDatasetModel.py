import os
import datetime
import json

from django.db import models
from django.core.validators import FileExtensionValidator
from django.dispatch import receiver

from accounts.models.CustomUserModel import CustomUserModel
from organisms.models.OrganismModel import OrganismModel

from ..storage import OverwriteStorage


def tenxfbcm_barcode_dataset_upload_directory(instance, filename):
    return 'tenxfbcm_datasets/{0}/{1}/{2}'.format(
        instance.tenxfbcm_dataset_upload_user.username,
        instance.tenxfbcm_dataset_folder,
        filename
    )


def tenxfbcm_feature_dataset_upload_directory(instance, filename):
    return 'tenxfbcm_datasets/{0}/{1}/{2}'.format(
        instance.tenxfbcm_dataset_upload_user.username,
        instance.tenxfbcm_dataset_folder,
        filename
    )


def tenxfbcm_matrix_dataset_upload_directory(instance, filename):
    return 'tenxfbcm_datasets/{0}/{1}/{2}'.format(
        instance.tenxfbcm_dataset_upload_user.username,
        instance.tenxfbcm_dataset_folder,
        filename
    )


class TenxfbcmDatasetModel(models.Model):
    tenxfbcm_barcode_dataset_file = models.FileField(
        upload_to=tenxfbcm_barcode_dataset_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['tsv.gz', 'tsv', 'gz'])]
    )
    tenxfbcm_feature_dataset_file = models.FileField(
        upload_to=tenxfbcm_feature_dataset_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['tsv.gz', 'tsv', 'gz'])]
    )
    tenxfbcm_matrix_dataset_file = models.FileField(
        upload_to=tenxfbcm_matrix_dataset_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['mtx.gz', 'mtx', 'gz'])]
    )
    tenxfbcm_dataset_name = models.CharField(max_length=200, null=False, blank=False)
    tenxfbcm_dataset_folder = models.CharField(max_length=100, null=False, blank=False)
    tenxfbcm_dataset_folder_path = models.CharField(max_length=200, null=True, blank=False)
    tenxfbcm_barcode_dataset_file_extension = models.CharField(max_length=50, null=False, blank=False)
    tenxfbcm_feature_dataset_file_extension = models.CharField(max_length=50, null=False, blank=False)
    tenxfbcm_matrix_dataset_file_extension = models.CharField(max_length=50, null=False, blank=False)
    tenxfbcm_dataset_organism = models.ForeignKey(OrganismModel, on_delete=models.CASCADE, null=False)
    tenxfbcm_dataset_public_flag = models.BooleanField(default=False)
    tenxfbcm_dataset_creation_timestamp = models.DateTimeField(auto_now_add=True)  # Automatically sets when created
    tenxfbcm_dataset_update_timestamp = models.DateTimeField(auto_now=True)  # Automatically updates on save
    tenxfbcm_dataset_upload_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)

    def __str__(self):
        output_string = json.dumps({
            "tenxfbcm_dataset_name": self.tenxfbcm_dataset_name
        })
        return output_string

    @property
    def get_tenxfbcm_dataset_name(self):
        return f"{self.tenxfbcm_dataset_name}"

    @property
    def get_tenxfbcm_dataset_public_flag(self):
        return f"{self.tenxfbcm_dataset_public_flag}"


@receiver(models.signals.post_save, sender=TenxfbcmDatasetModel)
def your_function_name(sender, instance, created, **kwargs):
    if created:
        if instance.tenxfbcm_dataset_folder:
            instance.tenxfbcm_dataset_folder_path = os.path.join(
                'tenxfbcm_datasets',
                instance.tenxfbcm_dataset_upload_user.username,
                instance.tenxfbcm_dataset_folder
            )
            instance.save()
    else:
        try:
            if instance.tenxfbcm_dataset_folder:
                if os.path.join('tenxfbcm_datasets', instance.tenxfbcm_dataset_upload_user.username, instance.tenxfbcm_dataset_folder) != instance.tenxfbcm_dataset_folder_path:
                    instance.tenxfbcm_dataset_folder_path = os.path.join(
                        'tenxfbcm_datasets',
                        instance.tenxfbcm_dataset_upload_user.username,
                        instance.tenxfbcm_dataset_folder
                    )
                    instance.save()
        except Exception as e:
            print(e)


@receiver(models.signals.post_delete, sender=TenxfbcmDatasetModel)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    tenxfbcm_barcode_dataset_file_path = instance.tenxfbcm_barcode_dataset_file.path
    tenxfbcm_feature_dataset_file_path = instance.tenxfbcm_feature_dataset_file.path
    tenxfbcm_matrix_dataset_file_path = instance.tenxfbcm_matrix_dataset_file.path
    if instance.tenxfbcm_barcode_dataset_file:
        instance.tenxfbcm_barcode_dataset_file.delete(save=False)
    if instance.tenxfbcm_feature_dataset_file:
        instance.tenxfbcm_feature_dataset_file.delete(save=False)
    if instance.tenxfbcm_matrix_dataset_file:
        instance.tenxfbcm_matrix_dataset_file.delete(save=False)
    tenxfbcm_folder_list = [
        os.path.dirname(os.path.dirname(tenxfbcm_barcode_dataset_file_path)),
        os.path.dirname(os.path.dirname(tenxfbcm_feature_dataset_file_path)),
        os.path.dirname(os.path.dirname(tenxfbcm_matrix_dataset_file_path))
    ]
    for i in range(0, len(tenxfbcm_folder_list)):
        try:
            if os.path.exists(tenxfbcm_folder_list[i]):
                os.rmdir(tenxfbcm_folder_list[i])
        except Exception as e:
            pass
