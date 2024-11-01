from django.contrib import admin

# Register your models here.
from .models.CustomUserModel import CustomUserModel

from .admins.CustomUserAdmin import CustomUserAdmin


# Register your models here.


admin.site.register(CustomUserModel, CustomUserAdmin)

