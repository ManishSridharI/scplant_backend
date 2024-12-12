import json

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.ModelModel import ModelModel

from ..serializers.ModelModelSerializer import ModelModelSerializer


@api_view(['GET'])
def ModelQuery(request):
    if request.method == 'GET':
        try:
            model_upload_user = request.user.id

            model = ModelModel.objects.filter(
                model_upload_user=model_upload_user
            )

            response_object = {
                "isModelQuery": True,
                "Model": ModelModelSerializer(model, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isModelQuery": False, "error": str(e)}, status=500)

    return Response({"isModelQuery": False, "error": "Invalid request method"}, status=405)
