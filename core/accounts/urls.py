from django.urls import re_path

from .api.TestAPI import Test


app_name = 'accounts'


urlpatterns = [
    re_path('api/test/', Test, name='api_test'),
]

