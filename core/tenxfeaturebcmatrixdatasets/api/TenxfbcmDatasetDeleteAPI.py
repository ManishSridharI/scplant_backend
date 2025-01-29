import json

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.TenxfbcmDatasetModel import TenxfbcmDatasetModel

from ..serializers.TenxfbcmDatasetModelSerializer import TenxfbcmDatasetModelSerializer


@api_view(['POST'])
def TenxfbcmDatasetDelete(request):
    if request.method == 'POST':
        try:
            tenxfbcm_dataset_id = request.data['tenxfbcm_dataset_id']

            tenxfbcm_dataset_instance = TenxfbcmDatasetModel.objects.filter(id=tenxfbcm_dataset_id)
            if tenxfbcm_dataset_instance:
                tenxfbcm_dataset_instance.delete()
                response_object = {"isTenxfbcmDatasetDelete": True}
                return Response(response_object)
            else:
                response_object = {
                    "isTenxfbcmDatasetDelete": False,
                    "error": "Data does not exists"
                }
                return Response(response_object, status=404)
        except Exception as e:
            return Response({"isTenxfbcmDatasetDelete": False, "error": str(e)}, status=500)

    return Response({"isTenxfbcmDatasetDelete": False, "error": "Invalid request method"}, status=405)
