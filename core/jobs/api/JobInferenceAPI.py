from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from celery.result import AsyncResult

from ..serializers.JobInferenceModelSerializer import JobInferenceModelSerializer
from ..serializers.JobInferenceArgumentModelSerializer import JobInferenceArgumentModelSerializer
from ..serializers.JobInferenceFileOutputModelSerializer import JobInferenceFileOutputModelSerializer

from ..tasks.InferenceTask import Inference

from scripts.models.ScriptModel import ScriptModel
from datasets.models.DatasetModel import DatasetModel
from predictors.models.PredictorModel import PredictorModel


@api_view(['POST'])
def JobInference(request):
    if request.method == 'POST':
        job_name = request.data['job_name']
        job_script = request.data['job_script']
        job_dataset = request.data['job_dataset']
        job_predictor = request.data['job_predictor']
        job_inference_gene_number = request.data['job_inference_gene_number']
        job_inference_log_filename = request.data['job_inference_log_filename']
        job_inference_prediction_filename = request.data['job_inference_prediction_filename']
        job_inference_stdout_filename = request.data['job_inference_stdout_filename']
        job_inference_stderr_filename = request.data['job_inference_stderr_filename']

        job_creation_user = request.user.id
        job_inference_file_creation_user = request.user.id
        job_inference_argument_creation_user = request.user.id

        try:
            script_instance = ScriptModel.objects.get(id=job_script)
        except Exception as e:
            return Response({"isJobInference": False, "error": "Script does not exists"}, status=404)

        try:
            dataset_instance = DatasetModel.objects.get(id=job_dataset)
        except Exception as e:
            return Response({"isJobInference": False, "error": "Dataset does not exists"}, status=404)

        try:
            predictor_instance = PredictorModel.objects.get(id=job_predictor)
        except Exception as e:
            return Response({"isJobInference": False, "error": "Predictor does not exists"}, status=404)

        try:
            job_inference_argument_model_serializer = JobInferenceArgumentModelSerializer(
                data={
                    "job_inference_gene_number": job_inference_gene_number,
                    "job_inference_argument_creation_user": job_inference_argument_creation_user
                }
            )

            job_inference_argument_model_serializer_instance = None
            if job_inference_argument_model_serializer.is_valid():
                try:
                    job_inference_argument_model_serializer_instance = job_inference_argument_model_serializer.save()
                except Exception as e:
                    return Response({"isJobInference": False, "error": str(e)}, status=405)
            else:
                return Response({"isJobInference": False, "error": str(job_inference_argument_model_serializer.errors)}, status=405)

            job_inference_file_output_model_serializer = JobInferenceFileOutputModelSerializer(
                data={
                    "job_inference_log_file": ContentFile("\n", name=str(job_inference_log_filename)+".txt"),
                    "job_inference_prediction_file": ContentFile("\n", name=str(job_inference_prediction_filename)+".csv"),
                    "job_inference_stdout_file": ContentFile("\n", name=str(job_inference_stdout_filename)+".txt"),
                    "job_inference_stderr_file": ContentFile("\n", name=str(job_inference_stderr_filename)+".txt"),
                    "job_inference_file_creation_user": job_inference_file_creation_user
                }
            )

            job_inference_file_output_model_serializer_instance = None
            if job_inference_file_output_model_serializer.is_valid():
                try:
                    job_inference_file_output_model_serializer_instance = job_inference_file_output_model_serializer.save()
                except Exception as e:
                    try:
                        if job_inference_argument_model_serializer_instance:
                            job_inference_argument_model_serializer_instance.delete()
                    except Exception as e:
                        return Response({"isJobInference": False, "error": str(e)}, status=405)
                    return Response({"isJobInference": False, "error": str(e)}, status=405)
            else:
                try:
                    if job_inference_argument_model_serializer_instance:
                        job_inference_argument_model_serializer_instance.delete()
                except Exception as e:
                    return Response({"isJobInference": False, "error": str(e)}, status=405)
                return Response({"isJobInference": False, "error": str(job_inference_file_output_model_serializer.errors)}, status=405)

            if job_inference_argument_model_serializer_instance and job_inference_file_output_model_serializer_instance:

                job_inference_argument = job_inference_argument_model_serializer_instance.id
                job_inference_file_output = job_inference_file_output_model_serializer_instance.id

                async_result_object = Inference.apply_async(
                    args=[
                        script_instance.script_file.path,
                        dataset_instance.dataset_file.path,
                        predictor_instance.predictor_file.path,
                        job_inference_argument_model_serializer_instance.job_inference_gene_number,
                        job_inference_file_output_model_serializer_instance.job_inference_log_file.path,
                        job_inference_file_output_model_serializer_instance.job_inference_prediction_file.path,
                        job_inference_file_output_model_serializer_instance.job_inference_stdout_file.path,
                        job_inference_file_output_model_serializer_instance.job_inference_stderr_file.path
                    ]
                )

                job_celery_task_id = async_result_object.id

                async_result_instance = AsyncResult(job_celery_task_id)

                job_inference_model_serializer = JobInferenceModelSerializer(
                    data={
                        "job_name": job_name,
                        "job_script": job_script,
                        "job_dataset": job_dataset,
                        "job_predictor": job_predictor,
                        "job_inference_argument": job_inference_argument,
                        "job_inference_file_output": job_inference_file_output,
                        "job_celery_task_id": job_celery_task_id,
                        "job_celery_task_status": async_result_instance.status,
                        "job_creation_user": job_creation_user
                    }
                )

                if job_inference_model_serializer.is_valid():
                    job_inference_model_serializer_instance = job_inference_model_serializer.save()

                    job_inference_argument_model_serializer = JobInferenceArgumentModelSerializer(
                        instance=job_inference_argument_model_serializer_instance,
                        data={
                            "job_celery_task_id": job_celery_task_id
                        },
                        partial=True
                    )

                    if job_inference_argument_model_serializer.is_valid():
                        try:
                            job_inference_argument_model_serializer_instance = job_inference_argument_model_serializer.save()
                        except Exception as e:
                            return Response({"isJobInference": False, "error": str(e)}, status=405)
                    else:
                        return Response({"isJobInference": False, "error": str(job_inference_argument_model_serializer.errors)}, status=405)

                    job_inference_file_output_model_serializer = JobInferenceFileOutputModelSerializer(
                        instance=job_inference_file_output_model_serializer_instance,
                        data={
                            "job_celery_task_id": job_celery_task_id
                        },
                        partial=True
                    )

                    if job_inference_file_output_model_serializer.is_valid():
                        try:
                            job_inference_file_output_model_serializer_instance = job_inference_file_output_model_serializer.save()
                        except Exception as e:
                            return Response({"isJobInference": False, "error": str(e)}, status=405)
                    else:
                        return Response({"isJobInference": False, "error": str(job_inference_file_output_model_serializer.errors)}, status=405)

                    return Response({"isJobInference": True}, status=201)
                else:
                    try:
                        if job_inference_argument_model_serializer_instance:
                            job_inference_argument_model_serializer_instance.delete()
                    except Exception as e:
                        return Response({"isJobInference": False, "error": str(e)}, status=405)
                    try:
                        if job_inference_file_output_model_serializer_instance:
                            job_inference_file_output_model_serializer_instance.delete()
                    except Exception as e:
                        return Response({"isJobInference": False, "error": str(e)}, status=405)
                    return Response({"isJobInference": False, "error": str(job_inference_model_serializer.errors)}, status=405)
            else:
                try:
                    if job_inference_argument_model_serializer_instance:
                        job_inference_argument_model_serializer_instance.delete()
                except Exception as e:
                    pass
                try:
                    if job_inference_file_output_model_serializer_instance:
                        job_inference_file_output_model_serializer_instance.delete()
                except Exception as e:
                    pass
                return Response({"isJobInference": False, "error": "Arguments and file outputs cannot be created"}, status=405)

        except Exception as e:
            return Response({"isJobInference": False, "error": str(e)}, status=405)

    return Response({"isJobInference": False, "error": "Invalid request method"}, status=405)
