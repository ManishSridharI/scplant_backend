import json

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.ModelModel import ModelModel

from ..serializers.ModelModelSerializer import ModelModelSerializer


@api_view(['POST'])
def ModelDelete(request):
    if request.method == 'POST':
        try:
            model_id = request.data['model_id']

            model_instance = ModelModel.objects.filter(id=model_id)
            if model_instance:
                model_instance.delete()
                response_object = {"isModelDelete": True}
                return Response(response_object)
            else:
                response_object = {
                    "isModelDelete": False,
                    "error": "Data does not exists"
                }
                return Response(response_object, status=404)
        except Exception as e:
            return Response({"isModelDelete": False, "error": str(e)}, status=500)

    return Response({"isModelDelete": False, "error": "Invalid request method"}, status=405)
