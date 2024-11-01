from django.db import models

from django.contrib.auth.models import AbstractUser


class CustomUserModel(AbstractUser):
    email = models.EmailField(null=False, blank=False, unique=True)
    first_name = models.CharField(max_length=200, null=False, blank=False)
    last_name = models.CharField(max_length=200, null=False, blank=True)
    organization = models.CharField(max_length=200, blank=True)
