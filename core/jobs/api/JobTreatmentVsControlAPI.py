from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from celery.result import AsyncResult

from ..serializers.JobTreatmentVsControlModelSerializer import JobTreatmentVsControlModelSerializer
from ..serializers.JobTreatmentVsControlFileOutputModelSerializer import JobTreatmentVsControlFileOutputModelSerializer

from ..tasks.TreatmentVsControlTask import TreatmentVsControl

from scripts.models.ScriptModel import ScriptModel
from datasets.models.DatasetModel import DatasetModel
from predictors.models.PredictorModel import PredictorModel


@api_view(['POST'])
def JobTreatmentVsControl(request):
    if request.method == 'POST':
        job_name = request.data['job_name']
        job_script = request.data['job_script']
        job_control_dataset = request.data['job_control_dataset']
        job_condition1_dataset = request.data['job_condition1_dataset']
        job_condition2_dataset = None
        try:
            job_condition2_dataset = request.data['job_condition2_dataset']
        except Exception as e:
            pass
        job_treatment_vs_control_stdout_filename = request.data['job_treatment_vs_control_stdout_filename']
        job_treatment_vs_control_stderr_filename = request.data['job_treatment_vs_control_stderr_filename']

        job_treatment_vs_control_condition1_marker_genes_filename = "control_vs_condition1_marker_genes"
        job_treatment_vs_control_condition1_top25_markers_filename = "control_vs_condition1_top25_markers"
        job_treatment_vs_control_condition1_top10_genes_dotplot_filename = "control_vs_condition1_top10_genes_dotplot"
        job_treatment_vs_control_condition2_marker_genes_filename = "control_vs_condition2_marker_genes"
        job_treatment_vs_control_condition2_top25_markers_filename = "control_vs_condition2_top25_markers"
        job_treatment_vs_control_condition2_top10_genes_dotplot_filename = "control_vs_condition2_top10_genes_dotplot"
        job_treatment_vs_control_control_vs_conditions_markers_filename = "control_vs_conditions_common_sig_markers"
        job_treatment_vs_control_conditions_vs_control_markers_filename = "conditions_vs_control_common_sig_markers"

        job_creation_user = request.user.id
        job_treatment_vs_control_file_creation_user = request.user.id

        try:
            script_instance = ScriptModel.objects.get(id=job_script)
        except Exception as e:
            return Response({"isJobTreatmentVsControl": False, "error": "Script does not exists"}, status=404)

        try:
            control_dataset_instance = DatasetModel.objects.get(id=job_control_dataset)
        except Exception as e:
            return Response({"isJobTreatmentVsControl": False, "error": "Control dataset does not exists"}, status=404)

        try:
            condition1_dataset_instance = DatasetModel.objects.get(id=job_condition1_dataset)
        except Exception as e:
            return Response({"isJobTreatmentVsControl": False, "error": "Condition1 dataset does not exists"}, status=404)

        condition2_dataset_instance = None
        if job_condition2_dataset:
            try:
                condition2_dataset_instance = DatasetModel.objects.get(id=job_condition2_dataset)
            except Exception as e:
                return Response({"isJobTreatmentVsControl": False, "error": "Condition2 dataset does not exists"}, status=404)

        try:
            job_treatment_vs_control_file_output_model_serializer = JobTreatmentVsControlFileOutputModelSerializer(
                data={
                    "job_treatment_vs_control_condition1_marker_genes_file": ContentFile("\n", name=str(job_treatment_vs_control_condition1_marker_genes_filename)+".csv"),
                    "job_treatment_vs_control_condition1_top25_markers_file": ContentFile("\n", name=str(job_treatment_vs_control_condition1_top25_markers_filename)+".txt"),
                    "job_treatment_vs_control_condition1_top10_genes_dotplot_file": ContentFile("\n", name=str(job_treatment_vs_control_condition1_top10_genes_dotplot_filename)+".pdf"),
                    "job_treatment_vs_control_condition2_marker_genes_file": ContentFile("\n", name=str(job_treatment_vs_control_condition2_marker_genes_filename)+".csv"),
                    "job_treatment_vs_control_condition2_top25_markers_file": ContentFile("\n", name=str(job_treatment_vs_control_condition2_top25_markers_filename)+".txt"),
                    "job_treatment_vs_control_condition2_top10_genes_dotplot_file": ContentFile("\n", name=str(job_treatment_vs_control_condition2_top10_genes_dotplot_filename)+".pdf"),
                    "job_treatment_vs_control_control_vs_conditions_markers_file": ContentFile("\n", name=str(job_treatment_vs_control_control_vs_conditions_markers_filename)+".txt"),
                    "job_treatment_vs_control_conditions_vs_control_markers_file": ContentFile("\n", name=str(job_treatment_vs_control_conditions_vs_control_markers_filename)+".txt"),
                    "job_treatment_vs_control_stdout_file": ContentFile("\n", name=str(job_treatment_vs_control_stdout_filename)+".txt"),
                    "job_treatment_vs_control_stderr_file": ContentFile("\n", name=str(job_treatment_vs_control_stderr_filename)+".txt"),
                    "job_treatment_vs_control_file_creation_user": job_treatment_vs_control_file_creation_user
                }
            )

            job_treatment_vs_control_file_output_model_serializer_instance = None
            if job_treatment_vs_control_file_output_model_serializer.is_valid():
                try:
                    job_treatment_vs_control_file_output_model_serializer_instance = job_treatment_vs_control_file_output_model_serializer.save()
                except Exception as e:
                    return Response({"isJobTreatmentVsControl": False, "error": str(e)}, status=405)
            else:
                return Response({"isJobTreatmentVsControl": False, "error": str(job_treatment_vs_control_file_output_model_serializer.errors)}, status=405)

            if job_treatment_vs_control_file_output_model_serializer_instance:

                job_treatment_vs_control_file_output = job_treatment_vs_control_file_output_model_serializer_instance.id

                treatment_vs_control_args = [
                    script_instance.script_file.path,
                    control_dataset_instance.dataset_file.path,
                    condition1_dataset_instance.dataset_file.path,
                    job_treatment_vs_control_file_output_model_serializer_instance.job_treatment_vs_control_condition1_marker_genes_file.path,
                    job_treatment_vs_control_file_output_model_serializer_instance.job_treatment_vs_control_condition1_top25_markers_file.path,
                    job_treatment_vs_control_file_output_model_serializer_instance.job_treatment_vs_control_condition1_top10_genes_dotplot_file.path,
                    job_treatment_vs_control_file_output_model_serializer_instance.job_treatment_vs_control_condition2_marker_genes_file.path,
                    job_treatment_vs_control_file_output_model_serializer_instance.job_treatment_vs_control_condition2_top25_markers_file.path,
                    job_treatment_vs_control_file_output_model_serializer_instance.job_treatment_vs_control_condition2_top10_genes_dotplot_file.path,
                    job_treatment_vs_control_file_output_model_serializer_instance.job_treatment_vs_control_control_vs_conditions_markers_file.path,
                    job_treatment_vs_control_file_output_model_serializer_instance.job_treatment_vs_control_conditions_vs_control_markers_file.path,
                    job_treatment_vs_control_file_output_model_serializer_instance.job_treatment_vs_control_stdout_file.path,
                    job_treatment_vs_control_file_output_model_serializer_instance.job_treatment_vs_control_stderr_file.path
                ]

                if condition2_dataset_instance:
                    treatment_vs_control_args += [condition2_dataset_instance.dataset_file.path]

                async_result_object = TreatmentVsControl.apply_async(
                    args=treatment_vs_control_args
                )

                job_celery_task_id = async_result_object.id

                async_result_instance = AsyncResult(job_celery_task_id)

                job_treatment_vs_control_model_serializer_data = {
                    "job_name": job_name,
                    "job_script": job_script,
                    "job_control_dataset": job_control_dataset,
                    "job_condition1_dataset": job_condition1_dataset,
                    "job_condition2_dataset": job_condition2_dataset,
                    "job_treatment_vs_control_file_output": job_treatment_vs_control_file_output,
                    "job_celery_task_id": job_celery_task_id,
                    "job_celery_task_status": async_result_instance.status,
                    "job_creation_user": job_creation_user
                }

                job_treatment_vs_control_model_serializer = JobTreatmentVsControlModelSerializer(
                    data=job_treatment_vs_control_model_serializer_data
                )

                if job_treatment_vs_control_model_serializer.is_valid():
                    job_treatment_vs_control_model_serializer_instance = job_treatment_vs_control_model_serializer.save()

                    job_treatment_vs_control_file_output_model_serializer = JobTreatmentVsControlFileOutputModelSerializer(
                        instance=job_treatment_vs_control_file_output_model_serializer_instance,
                        data={
                            "job_celery_task_id": job_celery_task_id
                        },
                        partial=True
                    )

                    if job_treatment_vs_control_file_output_model_serializer.is_valid():
                        try:
                            job_treatment_vs_control_file_output_model_serializer_instance = job_treatment_vs_control_file_output_model_serializer.save()
                        except Exception as e:
                            return Response({"isJobTreatmentVsControl": False, "error": str(e)}, status=405)
                    else:
                        return Response({"isJobTreatmentVsControl": False, "error": str(job_treatment_vs_control_file_output_model_serializer.errors)}, status=405)

                    return Response({"isJobTreatmentVsControl": True}, status=201)
                else:
                    try:
                        if job_treatment_vs_control_file_output_model_serializer_instance:
                            job_treatment_vs_control_file_output_model_serializer_instance.delete()
                    except Exception as e:
                        return Response({"isJobTreatmentVsControl": False, "error": str(e)}, status=405)
                    return Response({"isJobTreatmentVsControl": False, "error": str(job_treatment_vs_control_model_serializer.errors)}, status=405)
            else:
                try:
                    if job_treatment_vs_control_file_output_model_serializer_instance:
                        job_treatment_vs_control_file_output_model_serializer_instance.delete()
                except Exception as e:
                    pass
                return Response({"isJobTreatmentVsControl": False, "error": "File outputs cannot be created"}, status=405)

        except Exception as e:
            return Response({"isJobTreatmentVsControl": False, "error": str(e)}, status=405)

    return Response({"isJobTreatmentVsControl": False, "error": "Invalid request method"}, status=405)
