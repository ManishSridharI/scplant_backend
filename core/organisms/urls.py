from django.urls import re_path

from .api.OrganismGenerateAPI import OrganismGenerate
from .api.OrganismQueryAPI import OrganismQuery
from .api.OrganismDeleteAPI import OrganismDelete


app_name = 'organisms'


urlpatterns = [
    re_path('api/organism_generate/', OrganismGenerate, name='api_organism_generate'),
    re_path('api/organism_query/', OrganismQuery, name='api_organism_query'),
    re_path('api/organism_delete/', OrganismDelete, name='api_organism_delete')
]
