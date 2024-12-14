from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.PredictorModel import PredictorModel

from ..serializers.PredictorModelSerializer import PredictorModelSerializer


@api_view(['GET'])
def PredictorQuery(request):
    if request.method == 'GET':
        try:
            predictor_upload_user = request.user.id

            predictor = PredictorModel.objects.filter(
                predictor_upload_user=predictor_upload_user
            )

            response_object = {
                "isPredictorQuery": True,
                "Predictor": PredictorModelSerializer(predictor, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isPredictorQuery": False, "error": str(e)}, status=500)

    return Response({"isPredictorQuery": False, "error": "Invalid request method"}, status=405)


@api_view(['GET'])
def PredictorQueryPublic(request):
    if request.method == 'GET':
        try:
            predictor_upload_user = request.user.id

            predictor = PredictorModel.objects.filter(
                predictor_public_flag=True
            )

            response_object = {
                "isPredictorQueryPublic": True,
                "Predictor": PredictorModelSerializer(predictor, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isPredictorQueryPublic": False, "error": str(e)}, status=500)

    return Response({"isPredictorQueryPublic": False, "error": "Invalid request method"}, status=405)


@api_view(['GET'])
def PredictorQueryUploadedAndPublic(request):
    if request.method == 'GET':
        try:
            predictor_upload_user = request.user.id

            predictor = PredictorModel.objects.filter(
                Q(predictor_upload_user=predictor_upload_user) | Q(predictor_public_flag=True)
            )

            response_object = {
                "isPredictorQueryUploadedAndPublic": True,
                "Predictor": PredictorModelSerializer(predictor, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isPredictorQueryUploadedAndPublic": False, "error": str(e)}, status=500)

    return Response({"isPredictorQueryUploadedAndPublic": False, "error": "Invalid request method"}, status=405)
