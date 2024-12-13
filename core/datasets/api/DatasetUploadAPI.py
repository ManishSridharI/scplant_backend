from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.DatasetModel import DatasetModel

from ..serializers.DatasetModelSerializer import DatasetModelSerializer


@api_view(['POST'])
def DatasetUpload(request):
    if request.method == 'POST':
        dataset_file = request.FILES.get('dataset_file')
        dataset_name = request.POST.get('dataset_name')
        dataset_public_flag = request.POST.get('dataset_public_flag')
        dataset_upload_user = request.user.id

        serializer = DatasetModelSerializer(data={
            'dataset_name': dataset_name,
            'dataset_file': dataset_file,
            'dataset_public_flag': dataset_public_flag,
            'dataset_upload_user': dataset_upload_user
        })

        # Validate serializer and save if valid
        if serializer.is_valid():
            serializer.save()
            return Response({"isDatasetUploaded": True}, status=201)
        else:
            return Response({"isDatasetUploaded": False, "error": str(serializer.errors)}, status=405)

    return Response({"isDatasetUploaded": False, "error": "Invalid request method"}, status=405)
