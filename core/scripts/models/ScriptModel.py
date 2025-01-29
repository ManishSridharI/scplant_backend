import os
import datetime
import json

from django.db import models
from django.core.validators import FileExtensionValidator
from django.dispatch import receiver

from accounts.models.CustomUserModel import CustomUserModel


def script_upload_directory(instance, filename):
    return 'scripts/{0}/{1}/{2}'.format(instance.script_upload_user.username, str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%f")), filename)


class ScriptModel(models.Model):
    script_file = models.FileField(upload_to=script_upload_directory, validators=[FileExtensionValidator(allowed_extensions=['py', 'R'])])
    script_name = models.CharField(max_length=200, null=False, blank=False)
    script_file_extension = models.CharField(max_length=50, null=False, blank=False)
    script_public_flag = models.BooleanField(default=False)
    script_creation_timestamp = models.DateTimeField(auto_now_add=True)  # Automatically sets when created
    script_update_timestamp = models.DateTimeField(auto_now=True)  # Automatically updates on save
    script_upload_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)

    def __str__(self):
        output_string = json.dumps({
            "script_file": self.script_file,
            "script_name": self.script_name,
        })
        return output_string

    @property
    def get_script_name(self):
        return f"{self.script_name}"

    @property
    def get_script_path(self):
        return f"{self.script_path}"

    @property
    def get_script_public_flag(self):
        return f"{self.script_public_flag}"


@receiver(models.signals.post_delete, sender=ScriptModel)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    script_file_path = instance.script_file.path
    if instance.script_file:
        instance.script_file.delete(save=False)
    try:
        if os.path.exists(os.path.dirname(script_file_path)):
            os.rmdir(os.path.dirname(script_file_path))
    except Exception as e:
        pass
