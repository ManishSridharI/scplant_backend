from django.contrib.auth.admin import UserAdmin

from ..models.CustomUserModel import CustomUserModel


class CustomUserAdmin(UserAdmin):
	list_display = ['username', 'email', 'first_name', 'last_name', 'is_active', 'is_staff', 'is_superuser', 'organization']

	fieldsets = [
		(None, {
			'fields': ['username', 'password']
		}),
		('Personal info', {
			'fields': ['first_name', 'last_name', 'email', 'organization']
		}),
		('Permissions', {
			'fields': ['is_active', 'is_staff', 'is_superuser', 'groups', 'user_permissions']
		}),
		('Important dates', {
			'fields': ['last_login', 'date_joined']
		})
	]

	add_fieldsets = (
		(None, {
			'fields': ['username', 'password1', 'password2']
		}),
		('Personal info', {
			'fields': ['first_name', 'last_name', 'email', 'organization']
		}),
		('Permissions', {
			'fields': ['is_active', 'is_staff', 'is_superuser', 'groups', 'user_permissions']
		}),
		('Important dates', {
			'fields': ['last_login', 'date_joined']
		})
	)
