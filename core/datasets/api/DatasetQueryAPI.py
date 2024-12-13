from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.DatasetModel import DatasetModel

from ..serializers.DatasetModelSerializer import DatasetModelSerializer


@api_view(['GET'])
def DatasetQuery(request):
    if request.method == 'GET':
        try:
            dataset_upload_user = request.user.id

            dataset = DatasetModel.objects.filter(
                dataset_upload_user=dataset_upload_user
            )

            response_object = {
                "isDatasetQuery": True,
                "Dataset": DatasetModelSerializer(dataset, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isDatasetQuery": False, "error": str(e)}, status=500)

    return Response({"isDatasetQuery": False, "error": "Invalid request method"}, status=405)


@api_view(['GET'])
def DatasetQueryPublic(request):
    if request.method == 'GET':
        try:
            dataset_upload_user = request.user.id

            dataset = DatasetModel.objects.filter(
                dataset_public_flag=True
            )

            response_object = {
                "isDatasetQueryPublic": True,
                "Dataset": DatasetModelSerializer(dataset, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isDatasetQueryPublic": False, "error": str(e)}, status=500)

    return Response({"isDatasetQueryPublic": False, "error": "Invalid request method"}, status=405)


@api_view(['GET'])
def DatasetQueryUploadedAndPublic(request):
    if request.method == 'GET':
        try:
            dataset_upload_user = request.user.id

            dataset = DatasetModel.objects.filter(
                Q(dataset_upload_user=dataset_upload_user) | Q(dataset_public_flag=True)
            )

            response_object = {
                "isDatasetQueryUploadedAndPublic": True,
                "Dataset": DatasetModelSerializer(dataset, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isDatasetQueryUploadedAndPublic": False, "error": str(e)}, status=500)

    return Response({"isDatasetQueryUploadedAndPublic": False, "error": "Invalid request method"}, status=405)
