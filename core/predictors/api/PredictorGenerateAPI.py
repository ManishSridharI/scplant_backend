from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.PredictorModel import PredictorModel

from ..serializers.PredictorModelSerializer import PredictorModelSerializer


@api_view(['POST'])
def PredictorGenerate(request):
    if request.method == 'POST':
        predictor_name = request.data['predictor_name']
        predictor_filename = request.data['predictor_filename']
        predictor_organism = request.data['predictor_organism']
        predictor_public_flag = request.data['predictor_public_flag']
        predictor_upload_user = request.user.id

        predictor_model_serializer = PredictorModelSerializer(data={
            'predictor_name': predictor_name,
            'predictor_file': ContentFile("\n", name=str(predictor_filename)+".ckpt"),
            'predictor_organism': predictor_organism,
            'predictor_public_flag': predictor_public_flag,
            'predictor_upload_user': predictor_upload_user
        })

        # Validate serializer and save if valid
        if predictor_model_serializer.is_valid():
            predictor_model_serializer_instance = predictor_model_serializer.save()
            return Response({"isPredictorGenerate": True, "Predictor": predictor_model_serializer_instance.predictor_file.name}, status=201)
        else:
            return Response({"isPredictorGenerate": False, "error": str(predictor_model_serializer.errors)}, status=405)

    return Response({"isPredictorGenerate": False, "error": "Invalid request method"}, status=405)
