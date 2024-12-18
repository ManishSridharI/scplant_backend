import os
import datetime
import json

from django.db import models
from django.core.validators import FileExtensionValidator
from django.dispatch import receiver

from accounts.models.CustomUserModel import CustomUserModel


def dataset_upload_directory(instance, filename):
    return 'datasets/{0}/{1}/{2}'.format(instance.dataset_upload_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%f")), filename)


class DatasetModel(models.Model):
    dataset_file = models.FileField(upload_to=dataset_upload_directory, validators=[FileExtensionValidator(allowed_extensions=['h5ad'])])
    dataset_name = models.CharField(max_length=200, null=False, blank=False)
    dataset_public_flag = models.BooleanField(default=False)
    dataset_creation_timestamp = models.DateTimeField(auto_now_add=True)  # Automatically sets when created
    dataset_update_timestamp = models.DateTimeField(auto_now=True)  # Automatically updates on save
    dataset_upload_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)

    def __str__(self):
        output_string = json.dumps({
            "dataset_file": self.dataset_file,
            "dataset_name": self.dataset_name,
        })
        return output_string

    @property
    def get_dataset_name(self):
        return f"{self.dataset_name}"

    @property
    def get_dataset_path(self):
        return f"{self.dataset_path}"

    @property
    def get_dataset_public_flag(self):
        return f"{self.dataset_public_flag}"


@receiver(models.signals.post_delete, sender=DatasetModel)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    dataset_file_path = instance.dataset_file.path
    if instance.dataset_file:
        instance.dataset_file.delete(save=False)
    try:
        if os.path.exists(os.path.dirname(dataset_file_path)):
            os.rmdir(os.path.dirname(dataset_file_path))
    except Exception as e:
        pass
