import json

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
