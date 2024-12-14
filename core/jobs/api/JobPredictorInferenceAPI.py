from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from celery.result import AsyncResult

from ..models.JobModel import JobModel

from ..serializers.JobModelSerializer import JobModelSerializer

from ..tasks.Test import Add, WriteDate


@api_view(['POST'])
def JobPredictorInference(request):
    if request.method == 'POST':
        job_name = request.data['job_name']

        try:
            job_dataset = request.data['job_dataset']
        except Exception as e:
            job_dataset = None

        try:
            job_predictor = request.data['job_predictor']
        except Exception as e:
            job_predictor = None

        async_result_object = WriteDate.delay()

        job_celery_task = async_result_object.id

        job_creation_user = request.user.id

        serializer = JobModelSerializer(data={
            "job_name": job_name,
            "job_dataset": job_dataset,
            "job_predictor": job_predictor,
            "job_celery_task": job_celery_task,
            "job_creation_user": job_creation_user
        })

        # Validate serializer and save if valid
        if serializer.is_valid():
            serializer.save()
            return Response({"isJobPredictorInference": True}, status=201)
        else:
            return Response({"isJobPredictorInference": False, "error": str(serializer.errors)}, status=405)

    return Response({"isJobPredictorInference": False, "error": "Invalid request method"}, status=405)
