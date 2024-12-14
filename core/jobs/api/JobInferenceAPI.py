from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from celery.result import AsyncResult

from ..models.JobInferenceModel import JobInferenceModel

from ..serializers.JobInferenceModelSerializer import JobInferenceModelSerializer

from ..tasks.InferenceTask import Inference

from datasets.models.DatasetModel import DatasetModel
from predictors.models.PredictorModel import PredictorModel


@api_view(['POST'])
def JobInference(request):
    if request.method == 'POST':
        job_name = request.data['job_name']
        job_dataset = request.data['job_dataset']
        job_predictor = request.data['job_predictor']
        job_inference_gene_number = request.data['job_inference_gene_number']
        job_inference_log_filename = request.data['job_inference_log_filename']
        job_inference_prediction_filename = request.data['job_inference_prediction_filename']

        dataset_instance = DatasetModel.objects.get(id=job_dataset)
        predictor_instance = PredictorModel.objects.get(id=job_predictor)

        async_result_object = Inference.delay(
            job_inference_gene_number,
            dataset_instance.dataset_file.name,
            predictor_instance.predictor_file.name,
            str(job_inference_log_filename)+".txt",
            str(job_inference_prediction_filename)+".csv"
        )

        job_celery_task = async_result_object.id

        job_creation_user = request.user.id

        serializer = JobInferenceModelSerializer(data={
            "job_name": job_name,
            "job_dataset": job_dataset,
            "job_predictor": job_predictor,
            "job_celery_task": job_celery_task,
            "job_creation_user": job_creation_user,
            "job_inference_gene_number": job_inference_gene_number,
            "job_inference_log_file": ContentFile("\n", name=str(job_inference_log_filename)+".txt"),
            "job_inference_prediction_file": ContentFile("\n", name=str(job_inference_prediction_filename)+".csv")
        })

        # Validate serializer and save if valid
        if serializer.is_valid():
            serializer.save()
            return Response({"isJobInference": True}, status=201)
        else:
            return Response({"isJobInference": False, "error": str(serializer.errors)}, status=405)

    return Response({"isJobInference": False, "error": "Invalid request method"}, status=405)
