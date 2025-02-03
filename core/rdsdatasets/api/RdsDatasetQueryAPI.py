from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.RdsDatasetModel import RdsDatasetModel

from ..serializers.RdsDatasetModelSerializer import RdsDatasetModelSerializer


@api_view(['GET'])
def RdsDatasetQuery(request):
    if request.method == 'GET':
        try:
            rds_dataset_upload_user = request.user.id

            rds_dataset = RdsDatasetModel.objects.filter(
                rds_dataset_upload_user=rds_dataset_upload_user
            )

            response_object = {
                "isRdsDatasetQuery": True,
                "RdsDataset": RdsDatasetModelSerializer(rds_dataset, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isRdsDatasetQuery": False, "error": str(e)}, status=500)

    return Response({"isRdsDatasetQuery": False, "error": "Invalid request method"}, status=405)


@api_view(['GET'])
def RdsDatasetQueryPublic(request):
    if request.method == 'GET':
        try:
            rds_dataset_upload_user = request.user.id

            rds_dataset = RdsDatasetModel.objects.filter(
                rds_dataset_public_flag=True
            )

            response_object = {
                "isRdsDatasetQueryPublic": True,
                "RdsDataset": RdsDatasetModelSerializer(rds_dataset, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isRdsDatasetQueryPublic": False, "error": str(e)}, status=500)

    return Response({"isRdsDatasetQueryPublic": False, "error": "Invalid request method"}, status=405)


@api_view(['GET'])
def RdsDatasetQueryUploadedAndPublic(request):
    if request.method == 'GET':
        try:
            rds_dataset_upload_user = request.user.id

            rds_dataset = RdsDatasetModel.objects.filter(
                Q(rds_dataset_upload_user=rds_dataset_upload_user) | Q(rds_dataset_public_flag=True)
            )

            response_object = {
                "isRdsDatasetQueryUploadedAndPublic": True,
                "RdsDataset": RdsDatasetModelSerializer(rds_dataset, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isRdsDatasetQueryUploadedAndPublic": False, "error": str(e)}, status=500)

    return Response({"isRdsDatasetQueryUploadedAndPublic": False, "error": "Invalid request method"}, status=405)
