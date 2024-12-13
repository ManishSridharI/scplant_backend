from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.ModelModel import ModelModel

from ..serializers.ModelModelSerializer import ModelModelSerializer


@api_view(['POST'])
def ModelUpload(request):
    if request.method == 'POST':
        model_file = request.FILES.get('model_file')
        model_name = request.POST.get('model_name')
        model_public_flag = request.POST.get('model_public_flag')
        model_upload_user = request.user.id

        serializer = ModelModelSerializer(data={
            'model_name': model_name,
            'model_file': model_file,
            'model_public_flag': model_public_flag,
            'model_upload_user': model_upload_user
        })

        # Validate serializer and save if valid
        if serializer.is_valid():
            serializer.save()
            return Response({"isModelUploaded": True}, status=201)
        else:
            return Response({"isModelUploaded": False, "error": str(serializer.errors)}, status=405)

    return Response({"isModelUploaded": False, "error": "Invalid request method"}, status=405)
