from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.RdsDatasetModel import RdsDatasetModel

from ..serializers.RdsDatasetModelSerializer import RdsDatasetModelSerializer


@api_view(['POST'])
def RdsDatasetUpload(request):
    if request.method == 'POST':
        rds_dataset_file = request.FILES.get('rds_dataset_file')
        rds_dataset_name = request.POST.get('rds_dataset_name')
        rds_dataset_file_extension = request.POST.get('rds_dataset_file_extension')
        rds_dataset_organism = request.POST.get('rds_dataset_organism')
        rds_dataset_public_flag = request.POST.get('rds_dataset_public_flag')
        rds_dataset_upload_user = request.user.id

        rds_dataset_model_serializer = RdsDatasetModelSerializer(data={
            'rds_dataset_name': rds_dataset_name,
            'rds_dataset_file_extension': rds_dataset_file_extension,
            'rds_dataset_file': rds_dataset_file,
            'rds_dataset_organism': rds_dataset_organism,
            'rds_dataset_public_flag': rds_dataset_public_flag,
            'rds_dataset_upload_user': rds_dataset_upload_user
        })

        # Validate serializer and save if valid
        if rds_dataset_model_serializer.is_valid():
            rds_dataset_model_serializer_instance = rds_dataset_model_serializer.save()
            return Response({
                "isRdsDatasetUpload": True,
                "RdsDatasetUpload": {
                    "id": rds_dataset_model_serializer_instance.id,
                    "rds_dataset_name": rds_dataset_model_serializer_instance.rds_dataset_name,
                    "rds_dataset_file": rds_dataset_model_serializer_instance.rds_dataset_file.url,
                    "rds_dataset_file_extension": rds_dataset_model_serializer_instance.rds_dataset_file_extension,
                    "rds_dataset_organism": rds_dataset_model_serializer_instance.rds_dataset_organism.id,
                    "rds_dataset_public_flag": rds_dataset_model_serializer_instance.rds_dataset_public_flag,
                    "rds_dataset_upload_user": rds_dataset_model_serializer_instance.rds_dataset_upload_user.id
                }
            },
                status=201
            )
        else:
            return Response({"isRdsDatasetUpload": False, "error": str(rds_dataset_model_serializer.errors)}, status=405)

    return Response({"isRdsDatasetUpload": False, "error": "Invalid request method"}, status=405)
