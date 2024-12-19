from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from celery.result import AsyncResult

from ..serializers.JobCompareCellTypeDistModelSerializer import JobCompareCellTypeDistModelSerializer
from ..serializers.JobCompareCellTypeDistFileOutputModelSerializer import JobCompareCellTypeDistFileOutputModelSerializer

from ..tasks.CompareCellTypeDistTask import CompareCellTypeDist

from ..models.JobAnnotateAndPlotFileOutputModel import JobAnnotateAndPlotFileOutputModel

from scripts.models.ScriptModel import ScriptModel
from datasets.models.DatasetModel import DatasetModel
from predictors.models.PredictorModel import PredictorModel


@api_view(['POST'])
def JobCompareCellTypeDist(request):
    if request.method == 'POST':
        job_name = request.data['job_name']
        job_script = request.data['job_script']
        job_control_prediction_file = request.data['job_control_prediction_file']
        job_condition1_prediction_file = request.data['job_condition1_prediction_file']
        job_condition2_prediction_file = None
        try:
            job_condition2_prediction_file = request.data['job_condition2_prediction_file']
        except Exception as e:
            pass
        job_compare_cell_type_dist_stdout_filename = request.data['job_compare_cell_type_dist_stdout_filename']
        job_compare_cell_type_dist_stderr_filename = request.data['job_compare_cell_type_dist_stderr_filename']

        job_compare_cell_type_dist_output_filename = "compare_celltype_distributions"

        job_creation_user = request.user.id
        job_compare_cell_type_dist_file_creation_user = request.user.id

        try:
            script_instance = ScriptModel.objects.get(id=job_script)
        except Exception as e:
            return Response({"isJobCompareCellTypeDist": False, "error": "Script does not exists"}, status=404)

        try:
            control_prediction_file_instance = JobAnnotateAndPlotFileOutputModel.objects.get(id=job_control_prediction_file)
        except Exception as e:
            return Response({"isJobCompareCellTypeDist": False, "error": "Control dataset does not exists"}, status=404)

        try:
            condition1_prediction_file_instance = JobAnnotateAndPlotFileOutputModel.objects.get(id=job_condition1_prediction_file)
        except Exception as e:
            return Response({"isJobCompareCellTypeDist": False, "error": "Condition1 dataset does not exists"}, status=404)

        condition2_prediction_file_instance = None
        if job_condition2_prediction_file:
            try:
                condition2_prediction_file_instance = JobAnnotateAndPlotFileOutputModel.objects.get(id=job_condition2_prediction_file)
            except Exception as e:
                return Response({"isJobCompareCellTypeDist": False, "error": "Condition2 dataset does not exists"}, status=404)

        try:
            job_compare_cell_type_dist_file_output_model_serializer_instance = None
            job_compare_cell_type_dist_file_output_model_serializer = JobCompareCellTypeDistFileOutputModelSerializer(
                data={
                    "job_compare_cell_type_dist_output_file": ContentFile("\n", name=str(job_compare_cell_type_dist_output_filename)+".pdf"),
                    "job_compare_cell_type_dist_stdout_file": ContentFile("\n", name=str(job_compare_cell_type_dist_stdout_filename)+".txt"),
                    "job_compare_cell_type_dist_stderr_file": ContentFile("\n", name=str(job_compare_cell_type_dist_stderr_filename)+".txt"),
                    "job_compare_cell_type_dist_file_creation_user": job_compare_cell_type_dist_file_creation_user
                }
            )

            if job_compare_cell_type_dist_file_output_model_serializer.is_valid():
                try:
                    job_compare_cell_type_dist_file_output_model_serializer_instance = job_compare_cell_type_dist_file_output_model_serializer.save()
                except Exception as e:
                    return Response({"isJobCompareCellTypeDist": False, "error": str(e)}, status=405)
            else:
                return Response({"isJobCompareCellTypeDist": False, "error": str(job_compare_cell_type_dist_file_output_model_serializer.errors)}, status=405)

            if job_compare_cell_type_dist_file_output_model_serializer_instance:

                job_compare_cell_type_dist_file_output = job_compare_cell_type_dist_file_output_model_serializer_instance.id

                compare_cell_type_dist_args = [
                    script_instance.script_file.path,
                    control_prediction_file_instance.job_annotate_and_plot_prediction_file.path,
                    condition1_prediction_file_instance.job_annotate_and_plot_prediction_file.path,
                    job_compare_cell_type_dist_file_output_model_serializer_instance.job_compare_cell_type_dist_output_file.path,
                    job_compare_cell_type_dist_file_output_model_serializer_instance.job_compare_cell_type_dist_stdout_file.path,
                    job_compare_cell_type_dist_file_output_model_serializer_instance.job_compare_cell_type_dist_stderr_file.path
                ]
                if condition2_prediction_file_instance:
                    compare_cell_type_dist_args += [
                        condition2_prediction_file_instance.job_annotate_and_plot_prediction_file.path
                    ]

                async_result_object = CompareCellTypeDist.apply_async(
                    args=compare_cell_type_dist_args
                )

                job_celery_task_id = async_result_object.id

                async_result_instance = AsyncResult(job_celery_task_id)

                job_compare_cell_type_dist_model_serializer_data = {
                    "job_name": job_name,
                    "job_script": job_script,
                    "job_control_prediction_file": job_control_prediction_file,
                    "job_condition1_prediction_file": job_condition1_prediction_file,
                    "job_condition2_prediction_file": job_condition2_prediction_file,
                    "job_compare_cell_type_dist_file_output": job_compare_cell_type_dist_file_output,
                    "job_celery_task_id": job_celery_task_id,
                    "job_celery_task_status": async_result_instance.status,
                    "job_creation_user": job_creation_user
                }

                job_compare_cell_type_dist_model_serializer = JobCompareCellTypeDistModelSerializer(
                    data=job_compare_cell_type_dist_model_serializer_data
                )

                if job_compare_cell_type_dist_model_serializer.is_valid():
                    job_compare_cell_type_dist_model_serializer_instance = job_compare_cell_type_dist_model_serializer.save()

                    job_compare_cell_type_dist_file_output_model_serializer = JobCompareCellTypeDistFileOutputModelSerializer(
                        instance=job_compare_cell_type_dist_file_output_model_serializer_instance,
                        data={
                            "job_celery_task_id": job_celery_task_id
                        },
                        partial=True
                    )

                    if job_compare_cell_type_dist_file_output_model_serializer.is_valid():
                        try:
                            job_compare_cell_type_dist_file_output_model_serializer_instance = job_compare_cell_type_dist_file_output_model_serializer.save()
                        except Exception as e:
                            return Response({"isJobCompareCellTypeDist": False, "error": str(e)}, status=405)
                    else:
                        return Response({"isJobCompareCellTypeDist": False, "error": str(job_compare_cell_type_dist_file_output_model_serializer.errors)}, status=405)

                    return Response({"isJobCompareCellTypeDist": True}, status=201)
                else:
                    try:
                        if job_compare_cell_type_dist_file_output_model_serializer_instance:
                            job_compare_cell_type_dist_file_output_model_serializer_instance.delete()
                    except Exception as e:
                        return Response({"isJobCompareCellTypeDist": False, "error": str(e)}, status=405)
                    return Response({"isJobCompareCellTypeDist": False, "error": str(job_compare_cell_type_dist_model_serializer.errors)}, status=405)
            else:
                try:
                    if job_compare_cell_type_dist_file_output_model_serializer_instance:
                        job_compare_cell_type_dist_file_output_model_serializer_instance.delete()
                except Exception as e:
                    pass
                return Response({"isJobCompareCellTypeDist": False, "error": "File outputs cannot be created"}, status=405)

        except Exception as e:
            return Response({"isJobCompareCellTypeDist": False, "error": str(e)}, status=405)

    return Response({"isJobCompareCellTypeDist": False, "error": "Invalid request method"}, status=405)
