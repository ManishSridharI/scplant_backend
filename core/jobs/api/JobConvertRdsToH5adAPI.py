import os
import re

from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from celery.result import AsyncResult

from ..serializers.JobConvertRdsToH5adModelSerializer import JobConvertRdsToH5adModelSerializer
from ..serializers.JobConvertRdsToH5adFileOutputModelSerializer import JobConvertRdsToH5adFileOutputModelSerializer

from ..tasks.ConvertRdsToH5adTask import ConvertRdsToH5ad

from scripts.models.ScriptModel import ScriptModel
from rdsdatasets.models.RdsDatasetModel import RdsDatasetModel
from h5addatasets.models.H5adDatasetModel import H5adDatasetModel
from rdsdatasets.serializers.RdsDatasetModelSerializer import RdsDatasetModelSerializer
from h5addatasets.serializers.H5adDatasetModelSerializer import H5adDatasetModelSerializer
from predictors.models.PredictorModel import PredictorModel


@api_view(['POST'])
def JobConvertRdsToH5ad(request):
    if request.method == 'POST':
        job_name = request.data['job_name']
        job_script = request.data['job_script']
        job_rds_dataset = request.data['job_rds_dataset']
        job_convert_rds_to_h5ad_stdout_filename = request.data['job_convert_rds_to_h5ad_stdout_filename']
        job_convert_rds_to_h5ad_stderr_filename = request.data['job_convert_rds_to_h5ad_stderr_filename']

        job_creation_user = request.user.id
        job_convert_rds_to_h5ad_file_creation_user = request.user.id

        try:
            script_instance = ScriptModel.objects.get(id=job_script)
        except Exception as e:
            return Response({"isJobConvertRdsToH5ad": False, "error": "Script does not exists"}, status=404)

        try:
            rds_dataset_instance = RdsDatasetModel.objects.get(id=job_rds_dataset)
        except Exception as e:
            return Response({"isJobConvertRdsToH5ad": False, "error": "Dataset does not exists"}, status=404)

        try:
            job_convert_rds_to_h5ad_file_output_model_serializer_instance = None
            job_convert_rds_to_h5ad_file_output_model_serializer = JobConvertRdsToH5adFileOutputModelSerializer(
                data={
                    "job_organism": rds_dataset_instance.rds_dataset_organism.id,
                    "job_convert_rds_to_h5ad_output_file": ContentFile(
                        "\n",
                        name=str(
                            re.sub(
                                r'\.[^.]+$',
                                '',
                                os.path.basename(rds_dataset_instance.rds_dataset_file.path)
                            )
                        )+".h5ad"
                    ),
                    "job_convert_rds_to_h5ad_stdout_file": ContentFile("\n", name=str(job_convert_rds_to_h5ad_stdout_filename)+".txt"),
                    "job_convert_rds_to_h5ad_stderr_file": ContentFile("\n", name=str(job_convert_rds_to_h5ad_stderr_filename)+".txt"),
                    "job_convert_rds_to_h5ad_file_creation_user": job_convert_rds_to_h5ad_file_creation_user
                }
            )

            if job_convert_rds_to_h5ad_file_output_model_serializer.is_valid():
                try:
                    job_convert_rds_to_h5ad_file_output_model_serializer_instance = job_convert_rds_to_h5ad_file_output_model_serializer.save()
                except Exception as e:
                    return Response({"isJobConvertRdsToH5ad": False, "error": str(e)}, status=405)
            else:
                return Response({"isJobConvertRdsToH5ad": False, "error": str(job_convert_rds_to_h5ad_file_output_model_serializer.errors)}, status=405)


            h5ad_dataset_model_serializer_instance = None
            h5ad_dataset_model_serializer = H5adDatasetModelSerializer(
                data={
                    'h5ad_dataset_name': rds_dataset_instance.rds_dataset_name,
                    'h5ad_dataset_file_extension': "h5ad",
                    'h5ad_dataset_file': ContentFile("\n", name=re.sub(r'\.[^.]+$', '', str(os.path.basename(rds_dataset_instance.rds_dataset_file.path)))+"."+str("h5ad")),
                    'h5ad_dataset_organism': rds_dataset_instance.rds_dataset_organism.id,
                    'h5ad_dataset_public_flag': rds_dataset_instance.rds_dataset_public_flag,
                    'h5ad_dataset_upload_user': request.user.id
                }
            )

            if h5ad_dataset_model_serializer.is_valid():
                try:
                    h5ad_dataset_model_serializer_instance = h5ad_dataset_model_serializer.save()
                except Exception as e:
                    return Response({"isJobConvertRdsToH5ad": False, "error": str(e)}, status=405)
            else:
                try:
                    if job_convert_rds_to_h5ad_file_output_model_serializer_instance:
                        job_convert_rds_to_h5ad_file_output_model_serializer_instance.delete()
                except Exception as e:
                    pass
                return Response({"isJobConvertRdsToH5ad": False, "error": str(h5ad_dataset_model_serializer.errors)}, status=405)

            if h5ad_dataset_model_serializer_instance and job_convert_rds_to_h5ad_file_output_model_serializer_instance:

                job_convert_rds_to_h5ad_file_output = job_convert_rds_to_h5ad_file_output_model_serializer_instance.id

                convert_rds_to_h5ad_args = [
                    script_instance.script_file.path,
                    rds_dataset_instance.rds_dataset_file.path,
                    job_convert_rds_to_h5ad_file_output_model_serializer_instance.job_convert_rds_to_h5ad_output_file.path,
                    job_convert_rds_to_h5ad_file_output_model_serializer_instance.job_convert_rds_to_h5ad_stdout_file.path,
                    job_convert_rds_to_h5ad_file_output_model_serializer_instance.job_convert_rds_to_h5ad_stderr_file.path,
                    h5ad_dataset_model_serializer_instance.h5ad_dataset_file.path
                ]

                async_result_object = ConvertRdsToH5ad.apply_async(
                    args=convert_rds_to_h5ad_args
                )

                job_celery_task_id = async_result_object.id

                async_result_instance = AsyncResult(job_celery_task_id)

                job_convert_rds_to_h5ad_model_serializer_data = {
                    "job_name": job_name,
                    "job_script": job_script,
                    "job_rds_dataset": job_rds_dataset,
                    "job_organism": rds_dataset_instance.rds_dataset_organism.id,
                    "job_convert_rds_to_h5ad_file_output": job_convert_rds_to_h5ad_file_output,
                    "job_celery_task_id": job_celery_task_id,
                    "job_celery_task_status": async_result_instance.status,
                    "job_creation_user": job_creation_user
                }

                job_convert_rds_to_h5ad_model_serializer = JobConvertRdsToH5adModelSerializer(
                    data=job_convert_rds_to_h5ad_model_serializer_data
                )

                if job_convert_rds_to_h5ad_model_serializer.is_valid():
                    job_convert_rds_to_h5ad_model_serializer_instance = job_convert_rds_to_h5ad_model_serializer.save()

                    job_convert_rds_to_h5ad_file_output_model_serializer = JobConvertRdsToH5adFileOutputModelSerializer(
                        instance=job_convert_rds_to_h5ad_file_output_model_serializer_instance,
                        data={
                            "job_celery_task_id": job_celery_task_id
                        },
                        partial=True
                    )

                    if job_convert_rds_to_h5ad_file_output_model_serializer.is_valid():
                        try:
                            job_convert_rds_to_h5ad_file_output_model_serializer_instance = job_convert_rds_to_h5ad_file_output_model_serializer.save()
                        except Exception as e:
                            return Response({"isJobConvertRdsToH5ad": False, "error": str(e)}, status=405)
                    else:
                        return Response({"isJobConvertRdsToH5ad": False, "error": str(job_convert_rds_to_h5ad_file_output_model_serializer.errors)}, status=405)

                    return Response({"isJobConvertRdsToH5ad": True}, status=201)
                else:
                    try:
                        if job_convert_rds_to_h5ad_file_output_model_serializer_instance:
                            job_convert_rds_to_h5ad_file_output_model_serializer_instance.delete()
                    except Exception as e:
                        return Response({"isJobConvertRdsToH5ad": False, "error": str(e)}, status=405)
                    try:
                        if h5ad_dataset_model_serializer_instance:
                            h5ad_dataset_model_serializer_instance.delete()
                    except Exception as e:
                        return Response({"isJobConvertRdsToH5ad": False, "error": str(e)}, status=405)
                    return Response({"isJobConvertRdsToH5ad": False, "error": str(job_convert_rds_to_h5ad_model_serializer.errors)}, status=405)
            else:
                try:
                    if job_convert_rds_to_h5ad_file_output_model_serializer_instance:
                        job_convert_rds_to_h5ad_file_output_model_serializer_instance.delete()
                except Exception as e:
                    pass
                try:
                    if h5ad_dataset_model_serializer_instance:
                        h5ad_dataset_model_serializer_instance.delete()
                except Exception as e:
                    pass
                return Response({"isJobConvertRdsToH5ad": False, "error": "File outputs cannot be created"}, status=405)

        except Exception as e:
            return Response({"isJobConvertRdsToH5ad": False, "error": str(e)}, status=405)

    return Response({"isJobConvertRdsToH5ad": False, "error": "Invalid request method"}, status=405)
