import os
import datetime
import json

from django.db import models
from django.core.validators import FileExtensionValidator
from django.dispatch import receiver

from accounts.models.CustomUserModel import CustomUserModel
from organisms.models.OrganismModel import OrganismModel


def h5ad_dataset_upload_directory(instance, filename):
    return 'h5ad_datasets/{0}/{1}/{2}'.format(
        instance.h5ad_dataset_upload_user.username,
        str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%f")),
        filename
    )


class H5adDatasetModel(models.Model):
    h5ad_dataset_file = models.FileField(
        upload_to=h5ad_dataset_upload_directory,
        max_length=200,
        validators=[FileExtensionValidator(allowed_extensions=['h5ad'])]
    )
    h5ad_dataset_name = models.CharField(max_length=200, null=False, blank=False)
    h5ad_dataset_file_extension = models.CharField(max_length=50, null=False, blank=False)
    h5ad_dataset_organism = models.ForeignKey(OrganismModel, on_delete=models.CASCADE, null=False)
    h5ad_dataset_public_flag = models.BooleanField(default=False)
    h5ad_dataset_creation_timestamp = models.DateTimeField(auto_now_add=True)  # Automatically sets when created
    h5ad_dataset_update_timestamp = models.DateTimeField(auto_now=True)  # Automatically updates on save
    h5ad_dataset_upload_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)

    def __str__(self):
        output_string = json.dumps({
            "h5ad_dataset_name": self.h5ad_dataset_name
        })
        return output_string

    @property
    def get_h5ad_dataset_name(self):
        return f"{self.h5ad_dataset_name}"

    @property
    def get_h5ad_dataset_public_flag(self):
        return f"{self.h5ad_dataset_public_flag}"


@receiver(models.signals.post_delete, sender=H5adDatasetModel)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    h5ad_dataset_file_path = instance.h5ad_dataset_file.path
    if instance.h5ad_dataset_file:
        instance.h5ad_dataset_file.delete(save=False)
    try:
        if os.path.exists(os.path.dirname(h5ad_dataset_file_path)):
            os.rmdir(os.path.dirname(h5ad_dataset_file_path))
    except Exception as e:
        pass
