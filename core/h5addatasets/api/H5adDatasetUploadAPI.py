from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.H5adDatasetModel import H5adDatasetModel

from ..serializers.H5adDatasetModelSerializer import H5adDatasetModelSerializer


@api_view(['POST'])
def H5adDatasetUpload(request):
    if request.method == 'POST':
        h5ad_dataset_file = request.FILES.get('h5ad_dataset_file')
        h5ad_dataset_name = request.POST.get('h5ad_dataset_name')
        h5ad_dataset_file_extension = request.POST.get('h5ad_dataset_file_extension')
        h5ad_dataset_organism = request.POST.get('h5ad_dataset_organism')
        h5ad_dataset_public_flag = request.POST.get('h5ad_dataset_public_flag')
        h5ad_dataset_upload_user = request.user.id

        serializer = H5adDatasetModelSerializer(data={
            'h5ad_dataset_name': h5ad_dataset_name,
            'h5ad_dataset_file_extension': h5ad_dataset_file_extension,
            'h5ad_dataset_file': h5ad_dataset_file,
            'h5ad_dataset_organism': h5ad_dataset_organism,
            'h5ad_dataset_public_flag': h5ad_dataset_public_flag,
            'h5ad_dataset_upload_user': h5ad_dataset_upload_user
        })

        # Validate serializer and save if valid
        if serializer.is_valid():
            serializer.save()
            return Response({"isH5adDatasetUpload": True}, status=201)
        else:
            return Response({"isH5adDatasetUpload": False, "error": str(serializer.errors)}, status=405)

    return Response({"isH5adDatasetUpload": False, "error": "Invalid request method"}, status=405)
