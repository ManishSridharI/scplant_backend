from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.OrganismModel import OrganismModel

from ..serializers.OrganismModelSerializer import OrganismModelSerializer


@api_view(['POST'])
def OrganismGenerate(request):
    if request.method == 'POST':
        organism_name = request.data['organism_name']
        organism_creation_user = request.user.id

        organism_model_serializer = OrganismModelSerializer(data={
            'organism_name': organism_name,
            'organism_creation_user': organism_creation_user
        })

        # Validate serializer and save if valid
        if organism_model_serializer.is_valid():
            organism_model_serializer_instance = organism_model_serializer.save()
            return Response({"isOrganismGenerate": True, "Organism": organism_model_serializer_instance.organism_name}, status=201)
        else:
            return Response({"isOrganismGenerate": False, "error": str(organism_model_serializer.errors)}, status=405)

    return Response({"isOrganismGenerate": False, "error": "Invalid request method"}, status=405)
