from django.urls import re_path

from .api.ScriptGenerateAPI import ScriptGenerate
from .api.ScriptUploadAPI import ScriptUpload
from .api.ScriptQueryAPI import ScriptQuery, ScriptQueryPublic, ScriptQueryUploadedAndPublic
from .api.ScriptDeleteAPI import ScriptDelete


app_name = 'scripts'


urlpatterns = [
    re_path('api/script_generate/', ScriptGenerate, name='api_script_generate'),
    re_path('api/script_upload/', ScriptUpload, name='api_script_upload'),
    re_path('api/script_query/', ScriptQuery, name='api_script_query'),
    re_path('api/script_query_public/', ScriptQueryPublic, name='api_script_query_public'),
    re_path('api/script_query_uploaded_and_public/', ScriptQueryUploadedAndPublic, name='api_script_query_uploaded_and_public'),
    re_path('api/script_delete/', ScriptDelete, name='api_script_delete')
]
