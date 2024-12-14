from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.PredictorModel import PredictorModel

from ..serializers.PredictorModelSerializer import PredictorModelSerializer


@api_view(['POST'])
def PredictorUpload(request):
    if request.method == 'POST':
        predictor_file = request.FILES.get('predictor_file')
        predictor_name = request.POST.get('predictor_name')
        predictor_public_flag = request.POST.get('predictor_public_flag')
        predictor_upload_user = request.user.id

        serializer = PredictorModelSerializer(data={
            'predictor_name': predictor_name,
            'predictor_file': predictor_file,
            'predictor_public_flag': predictor_public_flag,
            'predictor_upload_user': predictor_upload_user
        })

        # Validate serializer and save if valid
        if serializer.is_valid():
            serializer.save()
            return Response({"isPredictorUpload": True}, status=201)
        else:
            return Response({"isPredictorUpload": False, "error": str(serializer.errors)}, status=405)

    return Response({"isPredictorUpload": False, "error": "Invalid request method"}, status=405)
