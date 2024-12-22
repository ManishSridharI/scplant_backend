from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.JobInferenceModel import JobInferenceModel

from ..serializers.JobInferenceModelSerializer import JobInferenceModelSerializer


@api_view(['GET'])
def JobInferenceQuery(request):
    if request.method == 'GET':
        try:
            job_creation_user = request.user.id

            job_inference_instance = JobInferenceModel.objects.filter(
                job_creation_user=job_creation_user
            )

            response_object = {
                "isJobInferenceQuery": True,
                "JobInference": JobInferenceModelSerializer(job_inference_instance, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isJobInferenceQuery": False, "error": str(e)}, status=500)

    return Response({"isJobInferenceQuery": False, "error": "Invalid request method"}, status=405)

