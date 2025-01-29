import os
import datetime
import json

from django.db import models
from django.core.validators import FileExtensionValidator
from django.dispatch import receiver

from accounts.models.CustomUserModel import CustomUserModel


class OrganismModel(models.Model):
    organism_name = models.CharField(max_length=200, null=False, blank=False)
    organism_creation_timestamp = models.DateTimeField(auto_now_add=True)  # Automatically sets when created
    organism_update_timestamp = models.DateTimeField(auto_now=True)  # Automatically updates on save
    organism_creation_user = models.ForeignKey(CustomUserModel, on_delete=models.CASCADE, null=False)

    def __str__(self):
        output_string = json.dumps({
            "organism_name": self.organism_name
        })
        return output_string

    @property
    def get_organism_name(self):
        return f"{self.organism_name}"
