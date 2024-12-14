import json

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.PredictorModel import PredictorModel

from ..serializers.PredictorModelSerializer import PredictorModelSerializer


@api_view(['POST'])
def PredictorDelete(request):
    if request.method == 'POST':
        try:
            predictor_id = request.data['predictor_id']

            predictor_instance = PredictorModel.objects.filter(id=predictor_id)
            if predictor_instance:
                predictor_instance.delete()
                response_object = {"isPredictorDelete": True}
                return Response(response_object)
            else:
                response_object = {
                    "isPredictorDelete": False,
                    "error": "Data does not exists"
                }
                return Response(response_object, status=404)
        except Exception as e:
            return Response({"isPredictorDelete": False, "error": str(e)}, status=500)

    return Response({"isPredictorDelete": False, "error": "Invalid request method"}, status=405)
