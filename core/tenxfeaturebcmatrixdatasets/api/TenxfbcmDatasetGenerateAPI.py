import datetime

from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.TenxfbcmDatasetModel import TenxfbcmDatasetModel

from ..serializers.TenxfbcmDatasetModelSerializer import TenxfbcmDatasetModelSerializer


@api_view(['POST'])
def TenxfbcmDatasetGenerate(request):
    if request.method == 'POST':
        tenxfbcm_dataset_name = request.data['tenxfbcm_dataset_name']
        tenxfbcm_dataset_folder = str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%f"))
        tenxfbcm_barcode_dataset_filename = request.data['tenxfbcm_barcode_dataset_filename']
        tenxfbcm_feature_dataset_filename = request.data['tenxfbcm_feature_dataset_filename']
        tenxfbcm_matrix_dataset_filename = request.data['tenxfbcm_matrix_dataset_filename']
        tenxfbcm_barcode_dataset_file_extension = request.data['tenxfbcm_barcode_dataset_file_extension']
        tenxfbcm_feature_dataset_file_extension = request.data['tenxfbcm_feature_dataset_file_extension']
        tenxfbcm_matrix_dataset_file_extension = request.data['tenxfbcm_matrix_dataset_file_extension']
        tenxfbcm_dataset_organism = request.data['tenxfbcm_dataset_organism']
        tenxfbcm_dataset_public_flag = request.data['tenxfbcm_dataset_public_flag']
        tenxfbcm_dataset_upload_user = request.user.id

        tenxfbcm_dataset_model_serializer = TenxfbcmDatasetModelSerializer(data={
            'tenxfbcm_dataset_name': tenxfbcm_dataset_name,
            'tenxfbcm_dataset_folder': tenxfbcm_dataset_folder,
            'tenxfbcm_barcode_dataset_file': ContentFile(
                "\n",
                name=str(tenxfbcm_barcode_dataset_filename) + "." + str(tenxfbcm_barcode_dataset_file_extension)
            ),
            'tenxfbcm_feature_dataset_file': ContentFile(
                "\n",
                name=str(tenxfbcm_feature_dataset_filename) + "." + str(tenxfbcm_feature_dataset_file_extension)
            ),
            'tenxfbcm_matrix_dataset_file': ContentFile(
                "\n",
                name=str(tenxfbcm_matrix_dataset_filename) + "." + str(tenxfbcm_matrix_dataset_file_extension)
            ),
            'tenxfbcm_barcode_dataset_file_extension': tenxfbcm_barcode_dataset_file_extension,
            'tenxfbcm_feature_dataset_file_extension': tenxfbcm_feature_dataset_file_extension,
            'tenxfbcm_matrix_dataset_file_extension': tenxfbcm_matrix_dataset_file_extension,
            'tenxfbcm_dataset_organism': tenxfbcm_dataset_organism,
            'tenxfbcm_dataset_public_flag': tenxfbcm_dataset_public_flag,
            'tenxfbcm_dataset_upload_user': tenxfbcm_dataset_upload_user
        })

        # Validate serializer and save if valid
        if tenxfbcm_dataset_model_serializer.is_valid():
            tenxfbcm_dataset_model_serializer_instance = tenxfbcm_dataset_model_serializer.save()
            return Response({
                "isTenxfbcmDatasetGenerate": True,
                "TenxfbcmDataset": tenxfbcm_dataset_model_serializer_instance.tenxfbcm_dataset_folder_path
            },
                status=201
            )
        else:
            return Response({
                "isTenxfbcmDatasetGenerate": False,
                "error": str(tenxfbcm_dataset_model_serializer.errors)
            },
                status=405
            )

    return Response({"isTenxfbcmDatasetGenerate": False, "error": "Invalid request method"}, status=405)
