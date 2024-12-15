from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.DatasetModel import DatasetModel

from ..serializers.DatasetModelSerializer import DatasetModelSerializer


@api_view(['POST'])
def DatasetGenerate(request):
    if request.method == 'POST':
        dataset_name = request.data['dataset_name']
        dataset_filename = request.data['dataset_filename']
        dataset_public_flag = request.data['dataset_public_flag']
        dataset_upload_user = request.user.id

        dataset_model_serializer = DatasetModelSerializer(data={
            'dataset_name': dataset_name,
            'dataset_file': ContentFile("\n", name=str(dataset_filename)+".h5ad"),
            'dataset_public_flag': dataset_public_flag,
            'dataset_upload_user': dataset_upload_user
        })

        # Validate serializer and save if valid
        if dataset_model_serializer.is_valid():
            dataset_model_serializer_instance = dataset_model_serializer.save()
            return Response({"isDatasetGenerate": True, "Dataset": dataset_model_serializer_instance.dataset_file.name}, status=201)
        else:
            return Response({"isDatasetGenerate": False, "error": str(dataset_model_serializer.errors)}, status=405)

    return Response({"isDatasetGenerate": False, "error": "Invalid request method"}, status=405)
