from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.H5adDatasetModel import H5adDatasetModel

from ..serializers.H5adDatasetModelSerializer import H5adDatasetModelSerializer


@api_view(['POST'])
def H5adDatasetGenerate(request):
    if request.method == 'POST':
        h5ad_dataset_name = request.data['h5ad_dataset_name']
        h5ad_dataset_file_extension = request.data['h5ad_dataset_file_extension']
        h5ad_dataset_filename = request.data['h5ad_dataset_filename']
        h5ad_dataset_organism = request.data['h5ad_dataset_organism']
        h5ad_dataset_public_flag = request.data['h5ad_dataset_public_flag']
        h5ad_dataset_upload_user = request.user.id

        h5ad_dataset_model_serializer = H5adDatasetModelSerializer(data={
            'h5ad_dataset_name': h5ad_dataset_name,
            'h5ad_dataset_file_extension': h5ad_dataset_file_extension,
            'h5ad_dataset_file': ContentFile("\n", name=str(h5ad_dataset_filename)+"."+str(h5ad_dataset_file_extension)),
            'h5ad_dataset_organism': h5ad_dataset_organism,
            'h5ad_dataset_public_flag': h5ad_dataset_public_flag,
            'h5ad_dataset_upload_user': h5ad_dataset_upload_user
        })

        # Validate serializer and save if valid
        if h5ad_dataset_model_serializer.is_valid():
            h5ad_dataset_model_serializer_instance = h5ad_dataset_model_serializer.save()
            return Response({"isH5adDatasetGenerate": True, "H5adDataset": h5ad_dataset_model_serializer_instance.h5ad_dataset_file.name}, status=201)
        else:
            return Response({"isH5adDatasetGenerate": False, "error": str(h5ad_dataset_model_serializer.errors)}, status=405)

    return Response({"isH5adDatasetGenerate": False, "error": "Invalid request method"}, status=405)
