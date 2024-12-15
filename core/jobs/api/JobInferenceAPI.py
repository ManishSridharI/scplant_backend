from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from celery.result import AsyncResult

from ..serializers.JobInferenceModelSerializer import JobInferenceModelSerializer
from ..serializers.JobInferenceFileOutputModelSerializer import JobInferenceFileOutputModelSerializer

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
        job_inference_stdout_filename = request.data['job_inference_stdout_filename']
        job_inference_stderr_filename = request.data['job_inference_stderr_filename']

        job_creation_user = request.user.id
        job_inference_file_creation_user = request.user.id

        try:
            dataset_instance = DatasetModel.objects.get(id=job_dataset)
        except Exception as e:
            return Response({"isJobInference": False, "error": "Data does not exists"}, status=404)

        try:
            predictor_instance = PredictorModel.objects.get(id=job_predictor)
        except Exception as e:
            return Response({"isJobInference": False, "error": "Data does not exists"}, status=404)

        job_inference_file_output_model_serializer = JobInferenceFileOutputModelSerializer(data={
            "job_inference_log_file": ContentFile("\n", name=str(job_inference_log_filename)+".txt"),
            "job_inference_prediction_file": ContentFile("\n", name=str(job_inference_prediction_filename)+".csv"),
            "job_inference_stdout_file": ContentFile("\n", name=str(job_inference_stdout_filename)+".txt"),
            "job_inference_stderr_file": ContentFile("\n", name=str(job_inference_stderr_filename)+".txt"),
            "job_inference_file_creation_user": job_inference_file_creation_user
        })

        # Validate serializer and save if valid
        if job_inference_file_output_model_serializer.is_valid():
            job_inference_file_output_model_serializer_instance = job_inference_file_output_model_serializer.save()

            job_inference_file_output = job_inference_file_output_model_serializer_instance.id

            async_result_object = Inference.delay(
                job_inference_gene_number,
                dataset_instance.dataset_file.path,
                predictor_instance.predictor_file.path,
                job_inference_file_output_model_serializer_instance.job_inference_log_file.path,
                job_inference_file_output_model_serializer_instance.job_inference_prediction_file.path,
                job_inference_file_output_model_serializer_instance.job_inference_stdout_file.path,
                job_inference_file_output_model_serializer_instance.job_inference_stderr_file.path
            )

            job_celery_task = async_result_object.id

            job_inference_model_serializer = JobInferenceModelSerializer(data={
                "job_name": job_name,
                "job_dataset": job_dataset,
                "job_predictor": job_predictor,
                "job_inference_gene_number": job_inference_gene_number,
                "job_inference_file_output": job_inference_file_output,
                "job_celery_task": job_celery_task,
                "job_creation_user": job_creation_user
            })

            if job_inference_model_serializer.is_valid():
                job_inference_model_serializer_instance = job_inference_model_serializer.save()
                return Response({"isJobInference": True}, status=201)
            else:
                return Response({"isJobInference": False, "error": str(job_inference_model_serializer.errors)}, status=405)
        else:
            return Response({"isJobInference": False, "error": str(job_inference_file_output_model_serializer.errors)}, status=405)

    return Response({"isJobInference": False, "error": "Invalid request method"}, status=405)
