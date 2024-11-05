from django.urls import re_path

from .api.TestAPI import Test
from .api.RegistrationAPI import Registration
from .api.PasswordChangeAPI import PasswordChange
from .api.LoginAPI import Login
from .api.LogoutAPI import Logout


app_name = 'accounts'


urlpatterns = [
    re_path('api/test/', Test, name='api_test'),
    re_path('api/registration', Registration, name='api_registration'),
    re_path('api/password_change/', PasswordChange, name='api_password_change'),
    re_path('api/login', Login, name='api_login'),
    re_path('api/logout', Logout, name='api_logout')
]

