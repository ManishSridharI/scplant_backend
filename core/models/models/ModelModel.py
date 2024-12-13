import json

from django.db import models
from django.core.validators import FileExtensionValidator

from accounts.models.CustomUserModel import CustomUserModel


class ModelModel(models.Model):
    model_file = models.FileField(upload_to='models/', validators=[FileExtensionValidator(allowed_extensions=['ckpt'])])
    model_name = models.CharField(max_length=200, null=False, blank=False)
    model_public_flag = models.BooleanField(default=False)
    model_creation_timestamp = models.DateTimeField(auto_now_add=True)  # Automatically sets when created
    model_update_timestamp = models.DateTimeField(auto_now=True)  # Automatically updates on save
    model_upload_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)

    def __str__(self):
        output_string = json.dumps({
            "model_file": self.model_file,
            "model_name": self.model_name,
            "model_public_flag": self.model_public_flag,
        })
        return output_string

    @property
    def get_model_name(self):
        return f"{self.model_name}"

    @property
    def get_model_path(self):
        return f"{self.model_path}"

    @property
    def get_model_public_flag(self):
        return f"{self.model_public_flag}"
