from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.OrganismModel import OrganismModel

from ..serializers.OrganismModelSerializer import OrganismModelSerializer


@api_view(['GET'])
def OrganismQuery(request):
    if request.method == 'GET':
        try:

            organism = OrganismModel.objects.all()

            response_object = {
                "isOrganismQuery": True,
                "Organism": OrganismModelSerializer(organism, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isOrganismQuery": False, "error": str(e)}, status=500)

    return Response({"isOrganismQuery": False, "error": "Invalid request method"}, status=405)
