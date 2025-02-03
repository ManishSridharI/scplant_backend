import os
import datetime

from django.conf import settings
from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from celery.result import AsyncResult

from ..serializers.JobTreatmentVsControlModelSerializer import JobTreatmentVsControlModelSerializer
from ..serializers.JobTreatmentVsControlFileOutputModelSerializer import JobTreatmentVsControlFileOutputModelSerializer

from ..tasks.TreatmentVsControlTask import TreatmentVsControl
from ..models.JobAnnotateAndPlotModel import JobAnnotateAndPlotModel
from ..models.JobAnnotateAndPlotFileOutputModel import JobAnnotateAndPlotFileOutputModel

from scripts.models.ScriptModel import ScriptModel
from h5addatasets.models.H5adDatasetModel import H5adDatasetModel
from tenxfeaturebcmatrixdatasets.models.TenxfbcmDatasetModel import TenxfbcmDatasetModel
from preddatasets.models.PredDatasetModel import PredDatasetModel
from predictors.models.PredictorModel import PredictorModel


@api_view(['POST'])
def JobTreatmentVsControl(request):
    if request.method == 'POST':
        job_name = request.data['job_name']
        job_script = request.data['job_script']

        job_control = request.data['job_control']
        job_condition1 = request.data['job_condition1']
        try:
            job_condition2 = request.data['job_condition2']
        except Exception as e:
            job_condition2 = None

        job_treatment_vs_control_data_type = None

        job_treatment_vs_control_folder = str(
            datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        )

        job_treatment_vs_control_stdout_filename = request.data['job_treatment_vs_control_stdout_filename']
        job_treatment_vs_control_stderr_filename = request.data['job_treatment_vs_control_stderr_filename']

        job_treatment_vs_control_comp_celltype_dist_filename = "compare_celltype_distributions"

        job_treatment_vs_control_control_vs_conditions_markers_filename = "control_vs_conditions_common_sig_markers"
        job_treatment_vs_control_conditions_vs_control_markers_filename = "conditions_vs_control_common_sig_markers"

        job_treatment_vs_control_ctrl_vs_cond1_cond1_all_DEGs_filename = "condition1_all_DEGs"
        job_treatment_vs_control_ctrl_vs_cond1_cond1_top10_DEGs_filename = "condition1_top10_DEGs"
        job_treatment_vs_control_ctrl_vs_cond1_cond1_top25_DEGs_filename = "condition1_top25_DEGs"
        job_treatment_vs_control_ctrl_vs_cond1_cond1_top5_DEGs_filename = "condition1_top5_DEGs"
        job_treatment_vs_control_ctrl_vs_cond1_ctrl_all_DEGs_filename = "control_all_DEGs"
        job_treatment_vs_control_ctrl_vs_cond1_ctrl_top10_DEGs_filename = "control_top10_DEGs"
        job_treatment_vs_control_ctrl_vs_cond1_ctrl_top25_DEGs_filename = "control_top25_DEGs"
        job_treatment_vs_control_ctrl_vs_cond1_ctrl_top5_DEGs_filename = "control_top5_DEGs"

        job_treatment_vs_control_ctrl_vs_cond2_cond2_all_DEGs_filename = "condition2_all_DEGs"
        job_treatment_vs_control_ctrl_vs_cond2_cond2_top10_DEGs_filename = "condition2_top10_DEGs"
        job_treatment_vs_control_ctrl_vs_cond2_cond2_top25_DEGs_filename = "condition2_top25_DEGs"
        job_treatment_vs_control_ctrl_vs_cond2_cond2_top5_DEGs_filename = "condition2_top5_DEGs"
        job_treatment_vs_control_ctrl_vs_cond2_ctrl_all_DEGs_filename = "control_all_DEGs"
        job_treatment_vs_control_ctrl_vs_cond2_ctrl_top10_DEGs_filename = "control_top10_DEGs"
        job_treatment_vs_control_ctrl_vs_cond2_ctrl_top25_DEGs_filename = "control_top25_DEGs"
        job_treatment_vs_control_ctrl_vs_cond2_ctrl_top5_DEGs_filename = "control_top5_DEGs"

        job_creation_user = request.user.id
        job_treatment_vs_control_file_creation_user = request.user.id

        control_job_organism = None
        condition1_job_organism = None
        condition2_job_organism = None
        control_data_path = None
        condition1_data_path = None
        condition2_data_path = None
        control_pred_file = None
        condition1_pred_file = None
        condition2_pred_file = None

        try:
            script_instance = ScriptModel.objects.get(id=job_script)
        except Exception as e:
            return Response({"isJobTreatmentVsControl": False, "error": "Script does not exists"}, status=404)

        try:
            control_dataset_instance = JobAnnotateAndPlotModel.objects.get(id=job_control)
            if control_dataset_instance:
                control_job_organism = control_dataset_instance.job_organism.id
                if control_dataset_instance.job_h5ad_dataset:
                    h5ad_dataset_instance = H5adDatasetModel.objects.get(
                        id=control_dataset_instance.job_h5ad_dataset.id
                    )
                    control_data_path = h5ad_dataset_instance.h5ad_dataset_file.path
                    job_treatment_vs_control_data_type = "h5ad"
                elif control_dataset_instance.job_tenxfbcm_dataset:
                    tenxfbcm_dataset_instance = TenxfbcmDatasetModel.objects.get(
                        id=control_dataset_instance.job_tenxfbcm_dataset.id
                    )
                    control_data_path = os.path.join(
                        settings.MEDIA_ROOT,
                        tenxfbcm_dataset_instance.tenxfbcm_dataset_folder_path
                    )
                    job_treatment_vs_control_data_type = "10x"
                else:
                    control_data_path = None
                if control_dataset_instance.job_annotate_and_plot_file_output:
                    job_annotate_and_plot_file_output_instance = JobAnnotateAndPlotFileOutputModel.objects.get(
                        id=control_dataset_instance.job_annotate_and_plot_file_output.id
                    )
                    control_pred_file = job_annotate_and_plot_file_output_instance.job_annotate_and_plot_prediction_file.path
                else:
                    control_pred_file = None
        except Exception as e:
            return Response({"isJobTreatmentVsControl": False, "error": str(e)}, status=404)

        try:
            condition1_dataset_instance = JobAnnotateAndPlotModel.objects.get(id=job_condition1)
            if condition1_dataset_instance:
                condition1_job_organism = control_dataset_instance.job_organism.id
                if condition1_dataset_instance.job_h5ad_dataset:
                    h5ad_dataset_instance = H5adDatasetModel.objects.get(
                        id=condition1_dataset_instance.job_h5ad_dataset.id
                    )
                    condition1_data_path = h5ad_dataset_instance.h5ad_dataset_file.path
                elif condition1_dataset_instance.job_tenxfbcm_dataset:
                    tenxfbcm_dataset_instance = TenxfbcmDatasetModel.objects.get(
                        id=condition1_dataset_instance.job_tenxfbcm_dataset.id
                    )
                    condition1_data_path = os.path.join(
                        settings.MEDIA_ROOT,
                        tenxfbcm_dataset_instance.tenxfbcm_dataset_folder_path
                    )
                else:
                    condition1_data_path = None
                if condition1_dataset_instance.job_annotate_and_plot_file_output:
                    job_annotate_and_plot_file_output_instance = JobAnnotateAndPlotFileOutputModel.objects.get(
                        id=condition1_dataset_instance.job_annotate_and_plot_file_output.id
                    )
                    condition1_pred_file = job_annotate_and_plot_file_output_instance.job_annotate_and_plot_prediction_file.path
                else:
                    condition1_pred_file = None
        except Exception as e:
            return Response({"isJobTreatmentVsControl": False, "error": str(e)}, status=404)

        if job_condition2:
            try:
                condition2_dataset_instance = JobAnnotateAndPlotModel.objects.get(id=job_condition2)
                if condition2_dataset_instance:
                    condition2_job_organism = condition2_dataset_instance.job_organism.id
                    if condition2_dataset_instance.job_h5ad_dataset:
                        h5ad_dataset_instance = H5adDatasetModel.objects.get(
                            id=condition2_dataset_instance.job_h5ad_dataset.id
                        )
                        condition2_data_path = h5ad_dataset_instance.h5ad_dataset_file.path
                    elif condition2_dataset_instance.job_tenxfbcm_dataset:
                        tenxfbcm_dataset_instance = TenxfbcmDatasetModel.objects.get(
                            id=condition2_dataset_instance.job_tenxfbcm_dataset.id
                        )
                        condition2_data_path = os.path.join(
                            settings.MEDIA_ROOT,
                            tenxfbcm_dataset_instance.tenxfbcm_dataset_folder_path
                        )
                    else:
                        condition2_data_path = None
                    if condition2_dataset_instance.job_annotate_and_plot_file_output:
                        job_annotate_and_plot_file_output_instance = JobAnnotateAndPlotFileOutputModel.objects.get(
                            id=condition2_dataset_instance.job_annotate_and_plot_file_output.id
                        )
                        condition2_pred_file = job_annotate_and_plot_file_output_instance.job_annotate_and_plot_prediction_file.path
                    else:
                        condition2_pred_file = None
            except Exception as e:
                return Response({"isJobTreatmentVsControl": False, "error": str(e)}, status=404)

        if control_job_organism != condition1_job_organism:
            return Response({"isJobTreatmentVsControl": False, "error": "Control and Condition1 organisms do not match"}, status=405)
        
        if job_condition2:
            if control_job_organism != condition2_job_organism:
                return Response({"isJobTreatmentVsControl": False, "error": "Control, Condition1, and Condition2 organisms do not match"}, status=405)

        job_organism = control_job_organism

        try:
            output_dict = {
                "job_organism": job_organism,
                "job_treatment_vs_control_folder": job_treatment_vs_control_folder,
                "job_treatment_vs_control_comp_celltype_dist_file": ContentFile("\n", name=str(job_treatment_vs_control_comp_celltype_dist_filename) + ".pdf"),
                "job_treatment_vs_control_ctrl_vs_cond1_cond1_all_DEGs_file": ContentFile("\n", name=str(job_treatment_vs_control_ctrl_vs_cond1_cond1_all_DEGs_filename) + ".xlsx"),
                "job_treatment_vs_control_ctrl_vs_cond1_cond1_top10_DEGs_file": ContentFile("\n", name=str(job_treatment_vs_control_ctrl_vs_cond1_cond1_top10_DEGs_filename) + ".xlsx"),
                "job_treatment_vs_control_ctrl_vs_cond1_cond1_top25_DEGs_file": ContentFile("\n", name=str(job_treatment_vs_control_ctrl_vs_cond1_cond1_top25_DEGs_filename) + ".xlsx"),
                "job_treatment_vs_control_ctrl_vs_cond1_cond1_top5_DEGs_file": ContentFile("\n", name=str(job_treatment_vs_control_ctrl_vs_cond1_cond1_top5_DEGs_filename) + ".xlsx"),
                "job_treatment_vs_control_ctrl_vs_cond1_ctrl_all_DEGs_file": ContentFile("\n", name=str(job_treatment_vs_control_ctrl_vs_cond1_ctrl_all_DEGs_filename) + ".xlsx"),
                "job_treatment_vs_control_ctrl_vs_cond1_ctrl_top10_DEGs_file": ContentFile("\n", name=str(job_treatment_vs_control_ctrl_vs_cond1_ctrl_top10_DEGs_filename) + ".xlsx"),
                "job_treatment_vs_control_ctrl_vs_cond1_ctrl_top25_DEGs_file": ContentFile("\n", name=str(job_treatment_vs_control_ctrl_vs_cond1_ctrl_top25_DEGs_filename) + ".xlsx"),
                "job_treatment_vs_control_ctrl_vs_cond1_ctrl_top5_DEGs_file": ContentFile("\n", name=str(job_treatment_vs_control_ctrl_vs_cond1_ctrl_top5_DEGs_filename) + ".xlsx"),
                "job_treatment_vs_control_stdout_file": ContentFile("\n", name=str(job_treatment_vs_control_stdout_filename) + ".txt"),
                "job_treatment_vs_control_stderr_file": ContentFile("\n", name=str(job_treatment_vs_control_stderr_filename) + ".txt"),
                "job_treatment_vs_control_file_creation_user": job_treatment_vs_control_file_creation_user
            }

            if job_condition2:
                output_dict["job_treatment_vs_control_control_vs_conditions_markers_file"] = ContentFile("\n", name=str(job_treatment_vs_control_control_vs_conditions_markers_filename) + ".xlsx")
                output_dict["job_treatment_vs_control_conditions_vs_control_markers_file"] = ContentFile("\n", name=str(job_treatment_vs_control_conditions_vs_control_markers_filename) + ".xlsx")
                output_dict["job_treatment_vs_control_ctrl_vs_cond2_cond2_all_DEGs_file"] = ContentFile("\n", name=str(job_treatment_vs_control_ctrl_vs_cond2_cond2_all_DEGs_filename) + ".xlsx")
                output_dict["job_treatment_vs_control_ctrl_vs_cond2_cond2_top10_DEGs_file"] = ContentFile("\n", name=str(job_treatment_vs_control_ctrl_vs_cond2_cond2_top10_DEGs_filename) + ".xlsx")
                output_dict["job_treatment_vs_control_ctrl_vs_cond2_cond2_top25_DEGs_file"] = ContentFile("\n", name=str(job_treatment_vs_control_ctrl_vs_cond2_cond2_top25_DEGs_filename) + ".xlsx")
                output_dict["job_treatment_vs_control_ctrl_vs_cond2_cond2_top5_DEGs_file"] = ContentFile("\n", name=str(job_treatment_vs_control_ctrl_vs_cond2_cond2_top5_DEGs_filename) + ".xlsx")
                output_dict["job_treatment_vs_control_ctrl_vs_cond2_ctrl_all_DEGs_file"] = ContentFile("\n", name=str(job_treatment_vs_control_ctrl_vs_cond2_ctrl_all_DEGs_filename) + ".xlsx")
                output_dict["job_treatment_vs_control_ctrl_vs_cond2_ctrl_top10_DEGs_file"] = ContentFile("\n", name=str(job_treatment_vs_control_ctrl_vs_cond2_ctrl_top10_DEGs_filename) + ".xlsx")
                output_dict["job_treatment_vs_control_ctrl_vs_cond2_ctrl_top25_DEGs_file"] = ContentFile("\n", name=str(job_treatment_vs_control_ctrl_vs_cond2_ctrl_top25_DEGs_filename) + ".xlsx")
                output_dict["job_treatment_vs_control_ctrl_vs_cond2_ctrl_top5_DEGs_file"] = ContentFile("\n", name=str(job_treatment_vs_control_ctrl_vs_cond2_ctrl_top5_DEGs_filename) + ".xlsx")

            job_treatment_vs_control_file_output_model_serializer_instance = None
            job_treatment_vs_control_file_output_model_serializer = JobTreatmentVsControlFileOutputModelSerializer(
                data=output_dict
            )

            if job_treatment_vs_control_file_output_model_serializer.is_valid():
                try:
                    job_treatment_vs_control_file_output_model_serializer_instance = job_treatment_vs_control_file_output_model_serializer.save()
                except Exception as e:
                    return Response({"isJobTreatmentVsControl": False, "error": str(e)}, status=405)
            else:
                return Response({
                    "isJobTreatmentVsControl": False,
                    "error": str(job_treatment_vs_control_file_output_model_serializer.errors)
                }, status=405)

            if job_treatment_vs_control_file_output_model_serializer_instance:

                job_treatment_vs_control_file_output = job_treatment_vs_control_file_output_model_serializer_instance.id

                treatment_vs_control_args = [
                    script_instance.script_file.path,
                    control_data_path,
                    condition1_data_path,
                    condition2_data_path,
                    control_pred_file,
                    condition1_pred_file,
                    condition2_pred_file,
                    os.path.join(
                        settings.MEDIA_ROOT,
                        job_treatment_vs_control_file_output_model_serializer_instance.job_treatment_vs_control_folder_path
                    ),
                    job_treatment_vs_control_file_output_model_serializer_instance.job_treatment_vs_control_stdout_file.path,
                    job_treatment_vs_control_file_output_model_serializer_instance.job_treatment_vs_control_stderr_file.path
                ]

                async_result_object = TreatmentVsControl.apply_async(
                    args=treatment_vs_control_args
                )

                job_celery_task_id = async_result_object.id

                async_result_instance = AsyncResult(job_celery_task_id)

                job_treatment_vs_control_model_serializer_data = {
                    "job_name": job_name,
                    "job_script": job_script,
                    "job_control": job_control,
                    "job_condition1": job_condition1,
                    "job_condition2": job_condition2,
                    "job_organism": job_organism,
                    "job_treatment_vs_control_data_type": job_treatment_vs_control_data_type,
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
                        return Response({
                            "isJobTreatmentVsControl": False,
                            "error": str(job_treatment_vs_control_file_output_model_serializer.errors)
                        },
                            status=405
                        )

                    return Response({"isJobTreatmentVsControl": True}, status=201)
                else:
                    try:
                        if job_treatment_vs_control_file_output_model_serializer_instance:
                            job_treatment_vs_control_file_output_model_serializer_instance.delete()
                    except Exception as e:
                        return Response({"isJobTreatmentVsControl": False, "error": str(e)}, status=405)
                    return Response({
                        "isJobTreatmentVsControl": False,
                        "error": str(job_treatment_vs_control_model_serializer.errors)
                    },
                        status=405
                    )
            else:
                try:
                    if job_treatment_vs_control_file_output_model_serializer_instance:
                        job_treatment_vs_control_file_output_model_serializer_instance.delete()
                except Exception as e:
                    pass
                return Response({
                    "isJobTreatmentVsControl": False,
                    "error": "File outputs cannot be created"
                },
                    status=405
                )

        except Exception as e:
            return Response({"isJobTreatmentVsControl": False, "error": str(e)}, status=405)

    return Response({"isJobTreatmentVsControl": False, "error": "Invalid request method"}, status=405)
