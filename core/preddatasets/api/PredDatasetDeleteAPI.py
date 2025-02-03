import json

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.PredDatasetModel import PredDatasetModel

from ..serializers.PredDatasetModelSerializer import PredDatasetModelSerializer


@api_view(['POST'])
def PredDatasetDelete(request):
    if request.method == 'POST':
        try:
            pred_dataset_id = request.data['pred_dataset_id']

            pred_dataset_instance = PredDatasetModel.objects.filter(id=pred_dataset_id)
            if pred_dataset_instance:
                pred_dataset_instance.delete()
                response_object = {"isPredDatasetDelete": True}
                return Response(response_object)
            else:
                response_object = {
                    "isPredDatasetDelete": False,
                    "error": "Data does not exists"
                }
                return Response(response_object, status=404)
        except Exception as e:
            return Response({"isPredDatasetDelete": False, "error": str(e)}, status=500)

    return Response({"isPredDatasetDelete": False, "error": "Invalid request method"}, status=405)
