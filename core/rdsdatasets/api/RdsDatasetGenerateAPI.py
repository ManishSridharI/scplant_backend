from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.RdsDatasetModel import RdsDatasetModel

from ..serializers.RdsDatasetModelSerializer import RdsDatasetModelSerializer


@api_view(['POST'])
def RdsDatasetGenerate(request):
    if request.method == 'POST':
        rds_dataset_name = request.data['rds_dataset_name']
        rds_dataset_file_extension = request.data['rds_dataset_file_extension']
        rds_dataset_filename = request.data['rds_dataset_filename']
        rds_dataset_organism = request.data['rds_dataset_organism']
        rds_dataset_public_flag = request.data['rds_dataset_public_flag']
        rds_dataset_upload_user = request.user.id

        rds_dataset_model_serializer = RdsDatasetModelSerializer(data={
            'rds_dataset_name': rds_dataset_name,
            'rds_dataset_file_extension': rds_dataset_file_extension,
            'rds_dataset_file': ContentFile("\n", name=str(rds_dataset_filename)+"."+str(rds_dataset_file_extension)),
            'rds_dataset_organism': rds_dataset_organism,
            'rds_dataset_public_flag': rds_dataset_public_flag,
            'rds_dataset_upload_user': rds_dataset_upload_user
        })

        # Validate serializer and save if valid
        if rds_dataset_model_serializer.is_valid():
            rds_dataset_model_serializer_instance = rds_dataset_model_serializer.save()
            return Response({"isRdsDatasetGenerate": True, "RdsDataset": rds_dataset_model_serializer_instance.rds_dataset_file.name}, status=201)
        else:
            return Response({"isRdsDatasetGenerate": False, "error": str(rds_dataset_model_serializer.errors)}, status=405)

    return Response({"isRdsDatasetGenerate": False, "error": "Invalid request method"}, status=405)
