import datetime

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.TenxfbcmDatasetModel import TenxfbcmDatasetModel

from ..serializers.TenxfbcmDatasetModelSerializer import TenxfbcmDatasetModelSerializer


@api_view(['POST'])
def TenxfbcmDatasetUpload(request):
    if request.method == 'POST':
        tenxfbcm_barcode_dataset_file = request.FILES.get('tenxfbcm_barcode_dataset_file')
        tenxfbcm_feature_dataset_file = request.FILES.get('tenxfbcm_feature_dataset_file')
        tenxfbcm_matrix_dataset_file = request.FILES.get('tenxfbcm_matrix_dataset_file')
        tenxfbcm_dataset_name = request.POST.get('tenxfbcm_dataset_name')
        tenxfbcm_dataset_folder = str(datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%f"))
        tenxfbcm_barcode_dataset_file_extension = request.POST.get('tenxfbcm_barcode_dataset_file_extension')
        tenxfbcm_feature_dataset_file_extension = request.POST.get('tenxfbcm_feature_dataset_file_extension')
        tenxfbcm_matrix_dataset_file_extension = request.POST.get('tenxfbcm_matrix_dataset_file_extension')
        tenxfbcm_dataset_organism = request.POST.get('tenxfbcm_dataset_organism')
        tenxfbcm_dataset_public_flag = request.POST.get('tenxfbcm_dataset_public_flag')
        tenxfbcm_dataset_upload_user = request.user.id

        serializer = TenxfbcmDatasetModelSerializer(data={
            'tenxfbcm_dataset_name': tenxfbcm_dataset_name,
            'tenxfbcm_dataset_folder': tenxfbcm_dataset_folder,
            'tenxfbcm_barcode_dataset_file': tenxfbcm_barcode_dataset_file,
            'tenxfbcm_feature_dataset_file': tenxfbcm_feature_dataset_file,
            'tenxfbcm_matrix_dataset_file': tenxfbcm_matrix_dataset_file,
            'tenxfbcm_barcode_dataset_file_extension': tenxfbcm_barcode_dataset_file_extension,
            'tenxfbcm_feature_dataset_file_extension': tenxfbcm_feature_dataset_file_extension,
            'tenxfbcm_matrix_dataset_file_extension': tenxfbcm_matrix_dataset_file_extension,
            'tenxfbcm_dataset_organism': tenxfbcm_dataset_organism,
            'tenxfbcm_dataset_public_flag': tenxfbcm_dataset_public_flag,
            'tenxfbcm_dataset_upload_user': tenxfbcm_dataset_upload_user
        })

        # Validate serializer and save if valid
        if serializer.is_valid():
            serializer.save()
            return Response({"isTenxfbcmDatasetUpload": True}, status=201)
        else:
            return Response({"isTenxfbcmDatasetUpload": False, "error": str(serializer.errors)}, status=405)

    return Response({"isTenxfbcmDatasetUpload": False, "error": "Invalid request method"}, status=405)
