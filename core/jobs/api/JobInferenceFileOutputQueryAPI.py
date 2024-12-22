from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.JobInferenceFileOutputModel import JobInferenceFileOutputModel

from ..serializers.JobInferenceFileOutputModelSerializer import JobInferenceFileOutputModelSerializer


@api_view(['GET'])
def JobInferenceFileOutputQuery(request):
    if request.method == 'GET':
        try:
            job_inference_file_creation_user = request.user.id

            job_inference_file_output_instance = JobInferenceFileOutputModel.objects.filter(
                job_inference_file_creation_user=job_inference_file_creation_user
            )

            response_object = {
                "isJobInferenceFileOutputQuery": True,
                "JobInferenceFileOutput": JobInferenceFileOutputModelSerializer(job_inference_file_output_instance, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isJobInferenceFileOutputQuery": False, "error": str(e)}, status=500)

    return Response({"isJobInferenceFileOutputQuery": False, "error": "Invalid request method"}, status=405)


@api_view(['POST'])
def JobInferenceFileOutputQueryByID(request):
    if request.method == 'POST':
        try:
            job_inference_file_output_id = request.data['job_inference_file_output_id']
            job_inference_file_creation_user = request.user.id

            job_inference_file_output_instance = JobInferenceFileOutputModel.objects.get(
                pk=job_inference_file_output_id,
                job_inference_file_creation_user=job_inference_file_creation_user
            )

            response_object = {
                "isJobInferenceFileOutputQueryByID": True,
                "JobInferenceFileOutput": JobInferenceFileOutputModelSerializer(job_inference_file_output_instance).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isJobInferenceFileOutputQueryByID": False, "error": str(e)}, status=500)

    return Response({"isJobInferenceFileOutputQueryByID": False, "error": "Invalid request method"}, status=405)
