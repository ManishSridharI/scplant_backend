import json

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.H5adDatasetModel import H5adDatasetModel

from ..serializers.H5adDatasetModelSerializer import H5adDatasetModelSerializer


@api_view(['POST'])
def H5adDatasetDelete(request):
    if request.method == 'POST':
        try:
            h5ad_dataset_id = request.data['h5ad_dataset_id']

            h5ad_dataset_instance = H5adDatasetModel.objects.filter(id=h5ad_dataset_id)
            if h5ad_dataset_instance:
                h5ad_dataset_instance.delete()
                response_object = {"isH5adDatasetDelete": True}
                return Response(response_object)
            else:
                response_object = {
                    "isH5adDatasetDelete": False,
                    "error": "Data does not exists"
                }
                return Response(response_object, status=404)
        except Exception as e:
            return Response({"isH5adDatasetDelete": False, "error": str(e)}, status=500)

    return Response({"isH5adDatasetDelete": False, "error": "Invalid request method"}, status=405)
