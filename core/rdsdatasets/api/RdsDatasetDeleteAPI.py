import json

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.RdsDatasetModel import RdsDatasetModel

from ..serializers.RdsDatasetModelSerializer import RdsDatasetModelSerializer


@api_view(['POST'])
def RdsDatasetDelete(request):
    if request.method == 'POST':
        try:
            rds_dataset_id = request.data['rds_dataset_id']

            rds_dataset_instance = RdsDatasetModel.objects.filter(id=rds_dataset_id)
            if rds_dataset_instance:
                rds_dataset_instance.delete()
                response_object = {"isRdsDatasetDelete": True}
                return Response(response_object)
            else:
                response_object = {
                    "isRdsDatasetDelete": False,
                    "error": "Data does not exists"
                }
                return Response(response_object, status=404)
        except Exception as e:
            return Response({"isRdsDatasetDelete": False, "error": str(e)}, status=500)

    return Response({"isRdsDatasetDelete": False, "error": "Invalid request method"}, status=405)
