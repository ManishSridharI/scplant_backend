from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.H5adDatasetModel import H5adDatasetModel

from ..serializers.H5adDatasetModelSerializer import H5adDatasetModelSerializer


@api_view(['GET'])
def H5adDatasetQuery(request):
    if request.method == 'GET':
        try:
            h5ad_dataset_upload_user = request.user.id

            h5ad_dataset = H5adDatasetModel.objects.filter(
                h5ad_dataset_upload_user=h5ad_dataset_upload_user
            )

            response_object = {
                "isH5adDatasetQuery": True,
                "H5adDataset": H5adDatasetModelSerializer(h5ad_dataset, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isH5adDatasetQuery": False, "error": str(e)}, status=500)

    return Response({"isH5adDatasetQuery": False, "error": "Invalid request method"}, status=405)


@api_view(['GET'])
def H5adDatasetQueryPublic(request):
    if request.method == 'GET':
        try:
            h5ad_dataset_upload_user = request.user.id

            h5ad_dataset = H5adDatasetModel.objects.filter(
                h5ad_dataset_public_flag=True
            )

            response_object = {
                "isH5adDatasetQueryPublic": True,
                "H5adDataset": H5adDatasetModelSerializer(h5ad_dataset, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isH5adDatasetQueryPublic": False, "error": str(e)}, status=500)

    return Response({"isH5adDatasetQueryPublic": False, "error": "Invalid request method"}, status=405)


@api_view(['GET'])
def H5adDatasetQueryUploadedAndPublic(request):
    if request.method == 'GET':
        try:
            h5ad_dataset_upload_user = request.user.id

            h5ad_dataset = H5adDatasetModel.objects.filter(
                Q(h5ad_dataset_upload_user=h5ad_dataset_upload_user) | Q(h5ad_dataset_public_flag=True)
            )

            response_object = {
                "isH5adDatasetQueryUploadedAndPublic": True,
                "H5adDataset": H5adDatasetModelSerializer(h5ad_dataset, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isH5adDatasetQueryUploadedAndPublic": False, "error": str(e)}, status=500)

    return Response({"isH5adDatasetQueryUploadedAndPublic": False, "error": "Invalid request method"}, status=405)
