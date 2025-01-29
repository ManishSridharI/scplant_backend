import json

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.OrganismModel import OrganismModel

from ..serializers.OrganismModelSerializer import OrganismModelSerializer


@api_view(['POST'])
def OrganismDelete(request):
    if request.method == 'POST':
        try:
            organism_id = request.data['organism_id']

            organism_instance = OrganismModel.objects.filter(id=organism_id)
            if organism_instance:
                organism_instance.delete()
                response_object = {"isOrganismDelete": True}
                return Response(response_object)
            else:
                response_object = {
                    "isOrganismDelete": False,
                    "error": "Data does not exists"
                }
                return Response(response_object, status=404)
        except Exception as e:
            return Response({"isOrganismDelete": False, "error": str(e)}, status=500)

    return Response({"isOrganismDelete": False, "error": "Invalid request method"}, status=405)
