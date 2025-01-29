from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.TenxfbcmDatasetModel import TenxfbcmDatasetModel

from ..serializers.TenxfbcmDatasetModelSerializer import TenxfbcmDatasetModelSerializer


@api_view(['GET'])
def TenxfbcmDatasetQuery(request):
    if request.method == 'GET':
        try:
            tenxfbcm_dataset_upload_user = request.user.id

            tenxfbcm_dataset = TenxfbcmDatasetModel.objects.filter(
                tenxfbcm_dataset_upload_user=tenxfbcm_dataset_upload_user
            )

            response_object = {
                "isTenxfbcmDatasetQuery": True,
                "TenxfbcmDataset": TenxfbcmDatasetModelSerializer(tenxfbcm_dataset, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isTenxfbcmDatasetQuery": False, "error": str(e)}, status=500)

    return Response({"isTenxfbcmDatasetQuery": False, "error": "Invalid request method"}, status=405)


@api_view(['GET'])
def TenxfbcmDatasetQueryPublic(request):
    if request.method == 'GET':
        try:
            tenxfbcm_dataset_upload_user = request.user.id

            tenxfbcm_dataset = TenxfbcmDatasetModel.objects.filter(
                tenxfbcm_dataset_public_flag=True
            )

            response_object = {
                "isTenxfbcmDatasetQueryPublic": True,
                "TenxfbcmDataset": TenxfbcmDatasetModelSerializer(tenxfbcm_dataset, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isTenxfbcmDatasetQueryPublic": False, "error": str(e)}, status=500)

    return Response({"isTenxfbcmDatasetQueryPublic": False, "error": "Invalid request method"}, status=405)


@api_view(['GET'])
def TenxfbcmDatasetQueryUploadedAndPublic(request):
    if request.method == 'GET':
        try:
            tenxfbcm_dataset_upload_user = request.user.id

            tenxfbcm_dataset = TenxfbcmDatasetModel.objects.filter(
                Q(tenxfbcm_dataset_upload_user=tenxfbcm_dataset_upload_user) | Q(tenxfbcm_dataset_public_flag=True)
            )

            response_object = {
                "isTenxfbcmDatasetQueryUploadedAndPublic": True,
                "TenxfbcmDataset": TenxfbcmDatasetModelSerializer(tenxfbcm_dataset, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isTenxfbcmDatasetQueryUploadedAndPublic": False, "error": str(e)}, status=500)

    return Response({"isTenxfbcmDatasetQueryUploadedAndPublic": False, "error": "Invalid request method"}, status=405)
