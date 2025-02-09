import os
import datetime
import json

from django.db import models
from django.core.validators import FileExtensionValidator
from django.dispatch import receiver

from accounts.models.CustomUserModel import CustomUserModel
from organisms.models.OrganismModel import OrganismModel


def predictor_upload_directory(instance, filename):
    return 'predictors/{0}/{1}/{2}'.format(instance.predictor_upload_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%f")), filename)


class PredictorModel(models.Model):
    predictor_file = models.FileField(upload_to=predictor_upload_directory, validators=[FileExtensionValidator(allowed_extensions=['ckpt'])])
    predictor_name = models.CharField(max_length=200, null=False, blank=False)
    predictor_organism = models.ForeignKey(OrganismModel, on_delete=models.CASCADE, null=False)
    predictor_public_flag = models.BooleanField(default=False)
    predictor_creation_timestamp = models.DateTimeField(auto_now_add=True)  # Automatically sets when created
    predictor_update_timestamp = models.DateTimeField(auto_now=True)  # Automatically updates on save
    predictor_upload_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)

    def __str__(self):
        output_string = json.dumps({
            "predictor_file": self.predictor_file,
            "predictor_name": self.predictor_name,
        })
        return output_string

    @property
    def get_predictor_name(self):
        return f"{self.predictor_name}"

    @property
    def get_predictor_path(self):
        return f"{self.predictor_path}"

    @property
    def get_predictor_public_flag(self):
        return f"{self.predictor_public_flag}"


@receiver(models.signals.post_delete, sender=PredictorModel)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    predictor_file_path = instance.predictor_file.path
    if instance.predictor_file:
        instance.predictor_file.delete(save=False)
    try:
        if os.path.exists(os.path.dirname(predictor_file_path)):
            os.rmdir(os.path.dirname(predictor_file_path))
    except Exception as e:
        pass
