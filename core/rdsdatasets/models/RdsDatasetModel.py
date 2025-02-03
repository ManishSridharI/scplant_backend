import os
import datetime
import json

from django.db import models
from django.core.validators import FileExtensionValidator
from django.dispatch import receiver

from accounts.models.CustomUserModel import CustomUserModel
from organisms.models.OrganismModel import OrganismModel

from ..storage import OverwriteStorage


def rds_dataset_upload_directory(instance, filename):
    return 'rds_datasets/{0}/{1}/{2}'.format(
        instance.rds_dataset_upload_user.username,
        str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%f")),
        filename
    )


class RdsDatasetModel(models.Model):
    rds_dataset_file = models.FileField(
        upload_to=rds_dataset_upload_directory,
        max_length=200,
        storage=OverwriteStorage(),
        validators=[FileExtensionValidator(allowed_extensions=['rds'])]
    )
    rds_dataset_name = models.CharField(max_length=200, null=False, blank=False)
    rds_dataset_file_extension = models.CharField(max_length=50, null=False, blank=False)
    rds_dataset_organism = models.ForeignKey(OrganismModel, on_delete=models.CASCADE, null=False)
    rds_dataset_public_flag = models.BooleanField(default=False)
    rds_dataset_creation_timestamp = models.DateTimeField(auto_now_add=True)  # Automatically sets when created
    rds_dataset_update_timestamp = models.DateTimeField(auto_now=True)  # Automatically updates on save
    rds_dataset_upload_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)

    def __str__(self):
        output_string = json.dumps({
            "rds_dataset_name": self.rds_dataset_name
        })
        return output_string

    @property
    def get_rds_dataset_name(self):
        return f"{self.rds_dataset_name}"

    @property
    def get_rds_dataset_public_flag(self):
        return f"{self.rds_dataset_public_flag}"


@receiver(models.signals.post_delete, sender=RdsDatasetModel)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    rds_dataset_file_path = instance.rds_dataset_file.path
    if instance.rds_dataset_file:
        instance.rds_dataset_file.delete(save=False)
    try:
        if os.path.exists(os.path.dirname(rds_dataset_file_path)):
            os.rmdir(os.path.dirname(rds_dataset_file_path))
    except Exception as e:
        pass
