from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.PredDatasetModel import PredDatasetModel

from ..serializers.PredDatasetModelSerializer import PredDatasetModelSerializer


@api_view(['POST'])
def PredDatasetGenerate(request):
    if request.method == 'POST':
        pred_dataset_name = request.data['pred_dataset_name']
        pred_dataset_file_extension = request.data['pred_dataset_file_extension']
        pred_dataset_filename = request.data['pred_dataset_filename']
        pred_dataset_organism = request.data['pred_dataset_organism']
        pred_dataset_public_flag = request.data['pred_dataset_public_flag']
        pred_dataset_upload_user = request.user.id

        pred_dataset_model_serializer = PredDatasetModelSerializer(data={
            'pred_dataset_name': pred_dataset_name,
            'pred_dataset_file_extension': pred_dataset_file_extension,
            'pred_dataset_file': ContentFile("\n", name=str(pred_dataset_filename)+"."+str(pred_dataset_file_extension)),
            'pred_dataset_organism': pred_dataset_organism,
            'pred_dataset_public_flag': pred_dataset_public_flag,
            'pred_dataset_upload_user': pred_dataset_upload_user
        })

        # Validate serializer and save if valid
        if pred_dataset_model_serializer.is_valid():
            pred_dataset_model_serializer_instance = pred_dataset_model_serializer.save()
            return Response({"isPredDatasetGenerate": True, "PredDataset": pred_dataset_model_serializer_instance.pred_dataset_file.name}, status=201)
        else:
            return Response({"isPredDatasetGenerate": False, "error": str(pred_dataset_model_serializer.errors)}, status=405)

    return Response({"isPredDatasetGenerate": False, "error": "Invalid request method"}, status=405)
