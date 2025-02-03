from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.PredDatasetModel import PredDatasetModel

from ..serializers.PredDatasetModelSerializer import PredDatasetModelSerializer


@api_view(['GET'])
def PredDatasetQuery(request):
    if request.method == 'GET':
        try:
            pred_dataset_upload_user = request.user.id

            pred_dataset = PredDatasetModel.objects.filter(
                pred_dataset_upload_user=pred_dataset_upload_user
            )

            response_object = {
                "isPredDatasetQuery": True,
                "PredDataset": PredDatasetModelSerializer(pred_dataset, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isPredDatasetQuery": False, "error": str(e)}, status=500)

    return Response({"isPredDatasetQuery": False, "error": "Invalid request method"}, status=405)


@api_view(['GET'])
def PredDatasetQueryPublic(request):
    if request.method == 'GET':
        try:
            pred_dataset_upload_user = request.user.id

            pred_dataset = PredDatasetModel.objects.filter(
                pred_dataset_public_flag=True
            )

            response_object = {
                "isPredDatasetQueryPublic": True,
                "PredDataset": PredDatasetModelSerializer(pred_dataset, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isPredDatasetQueryPublic": False, "error": str(e)}, status=500)

    return Response({"isPredDatasetQueryPublic": False, "error": "Invalid request method"}, status=405)


@api_view(['GET'])
def PredDatasetQueryUploadedAndPublic(request):
    if request.method == 'GET':
        try:
            pred_dataset_upload_user = request.user.id

            pred_dataset = PredDatasetModel.objects.filter(
                Q(pred_dataset_upload_user=pred_dataset_upload_user) | Q(pred_dataset_public_flag=True)
            )

            response_object = {
                "isPredDatasetQueryUploadedAndPublic": True,
                "PredDataset": PredDatasetModelSerializer(pred_dataset, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isPredDatasetQueryUploadedAndPublic": False, "error": str(e)}, status=500)

    return Response({"isPredDatasetQueryUploadedAndPublic": False, "error": "Invalid request method"}, status=405)
