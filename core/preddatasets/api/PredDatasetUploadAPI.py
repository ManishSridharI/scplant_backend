from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.PredDatasetModel import PredDatasetModel

from ..serializers.PredDatasetModelSerializer import PredDatasetModelSerializer


@api_view(['POST'])
def PredDatasetUpload(request):
    if request.method == 'POST':
        pred_dataset_file = request.FILES.get('pred_dataset_file')
        pred_dataset_name = request.POST.get('pred_dataset_name')
        pred_dataset_file_extension = request.POST.get('pred_dataset_file_extension')
        pred_dataset_organism = request.POST.get('pred_dataset_organism')
        pred_dataset_public_flag = request.POST.get('pred_dataset_public_flag')
        pred_dataset_upload_user = request.user.id

        serializer = PredDatasetModelSerializer(data={
            'pred_dataset_name': pred_dataset_name,
            'pred_dataset_file_extension': pred_dataset_file_extension,
            'pred_dataset_file': pred_dataset_file,
            'pred_dataset_organism': pred_dataset_organism,
            'pred_dataset_public_flag': pred_dataset_public_flag,
            'pred_dataset_upload_user': pred_dataset_upload_user
        })

        # Validate serializer and save if valid
        if serializer.is_valid():
            serializer.save()
            return Response({"isPredDatasetUpload": True}, status=201)
        else:
            return Response({"isPredDatasetUpload": False, "error": str(serializer.errors)}, status=405)

    return Response({"isPredDatasetUpload": False, "error": "Invalid request method"}, status=405)
