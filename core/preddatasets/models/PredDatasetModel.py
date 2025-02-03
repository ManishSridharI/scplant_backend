import os
import datetime
import json

from django.db import models
from django.core.validators import FileExtensionValidator
from django.dispatch import receiver

from accounts.models.CustomUserModel import CustomUserModel
from organisms.models.OrganismModel import OrganismModel


def pred_dataset_upload_directory(instance, filename):
    return 'pred_datasets/{0}/{1}/{2}'.format(
        instance.pred_dataset_upload_user.username,
        str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%f")),
        filename
    )


class PredDatasetModel(models.Model):
    pred_dataset_file = models.FileField(
        upload_to=pred_dataset_upload_directory,
        max_length=200,
        validators=[FileExtensionValidator(allowed_extensions=['csv'])]
    )
    pred_dataset_name = models.CharField(max_length=200, null=False, blank=False)
    pred_dataset_file_extension = models.CharField(max_length=50, null=False, blank=False)
    pred_dataset_organism = models.ForeignKey(OrganismModel, on_delete=models.CASCADE, null=False)
    pred_dataset_public_flag = models.BooleanField(default=False)
    pred_dataset_creation_timestamp = models.DateTimeField(auto_now_add=True)  # Automatically sets when created
    pred_dataset_update_timestamp = models.DateTimeField(auto_now=True)  # Automatically updates on save
    pred_dataset_upload_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)

    def __str__(self):
        output_string = json.dumps({
            "pred_dataset_name": self.pred_dataset_name
        })
        return output_string

    @property
    def get_pred_dataset_name(self):
        return f"{self.pred_dataset_name}"

    @property
    def get_pred_dataset_public_flag(self):
        return f"{self.pred_dataset_public_flag}"


@receiver(models.signals.post_delete, sender=PredDatasetModel)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    pred_dataset_file_path = instance.pred_dataset_file.path
    if instance.pred_dataset_file:
        instance.pred_dataset_file.delete(save=False)
    try:
        if os.path.exists(os.path.dirname(pred_dataset_file_path)):
            os.rmdir(os.path.dirname(pred_dataset_file_path))
    except Exception as e:
        pass
