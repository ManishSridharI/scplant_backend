import json

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.DatasetModel import DatasetModel

from ..serializers.DatasetModelSerializer import DatasetModelSerializer


@api_view(['POST'])
def DatasetDelete(request):
    if request.method == 'POST':
        try:
            dataset_id = request.data['dataset_id']

            dataset_instance = DatasetModel.objects.filter(id=dataset_id)
            if dataset_instance:
                dataset_instance.delete()
                response_object = {"isDatasetDelete": True}
                return Response(response_object)
            else:
                response_object = {
                    "isDatasetDelete": False,
                    "error": "Data does not exists"
                }
                return Response(response_object, status=404)
        except Exception as e:
            return Response({"isDatasetDelete": False, "error": str(e)}, status=500)

    return Response({"isDatasetDelete": False, "error": "Invalid request method"}, status=405)
