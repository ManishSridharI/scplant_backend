import os
import datetime

from django.conf import settings
from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from celery.result import AsyncResult

from ..serializers.JobAnnotateAndPlotModelSerializer import JobAnnotateAndPlotModelSerializer
from ..serializers.JobAnnotateAndPlotFileOutputModelSerializer import JobAnnotateAndPlotFileOutputModelSerializer

from ..tasks.AnnotateAndPlotTask import AnnotateAndPlot

from scripts.models.ScriptModel import ScriptModel
from h5addatasets.models.H5adDatasetModel import H5adDatasetModel
from tenxfeaturebcmatrixdatasets.models.TenxfbcmDatasetModel import TenxfbcmDatasetModel
from preddatasets.models.PredDatasetModel import PredDatasetModel
from predictors.models.PredictorModel import PredictorModel
from preddatasets.serializers.PredDatasetModelSerializer import PredDatasetModelSerializer


@api_view(['POST'])
def JobAnnotateAndPlot(request):
    if request.method == 'POST':
        job_name = request.data['job_name']
        job_script = request.data['job_script']

        try:
            job_h5ad_dataset = request.data['job_h5ad_dataset']
        except Exception as e:
            job_h5ad_dataset = None
        try:
            job_tenxfbcm_dataset = request.data['job_tenxfbcm_dataset']
        except Exception as e:
            job_tenxfbcm_dataset = None
        try:
            job_pred_dataset = request.data['job_pred_dataset']
        except Exception as e:
            job_pred_dataset = None

        job_predictor = request.data['job_predictor']

        job_annotate_and_plot_data_type = None

        job_annotate_and_plot_folder = str(
            datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        )

        job_annotate_and_plot_stdout_filename = request.data['job_annotate_and_plot_stdout_filename']
        job_annotate_and_plot_stderr_filename = request.data['job_annotate_and_plot_stderr_filename']

        job_annotate_and_plot_output_with_celltype_filename = "output_with_celltype"
        job_annotate_and_plot_prediction_filename = "prediction"
        job_annotate_and_plot_log_filename = "log"
        job_annotate_and_plot_stats_csv_filename = "stats"
        job_annotate_and_plot_stats_pdf_filename = "stats"
        job_annotate_and_plot_marker_genes_filename = "marker_genes"
        job_annotate_and_plot_annotate_tsne_filename = "annotate_tsne"
        job_annotate_and_plot_annotate_umap_filename = "annotate_umap"
        job_annotate_and_plot_top3_genes_dotplot_filename = "top3_genes_dotplot"
        job_annotate_and_plot_all_markers_filename = "all_marker_genes"
        job_annotate_and_plot_top5_markers_filename = "top5_marker_genes"
        job_annotate_and_plot_top10_markers_filename = "top10_marker_genes"
        job_annotate_and_plot_top25_markers_filename = "top25_marker_genes"

        job_creation_user = request.user.id
        job_annotate_and_plot_file_creation_user = request.user.id
        pred_dataset_upload_user = request.user.id

        try:
            script_instance = ScriptModel.objects.get(id=job_script)
        except Exception as e:
            return Response({"isJobAnnotateAndPlot": False, "error": "Script does not exists"}, status=404)

        try:
            h5ad_dataset_instance = H5adDatasetModel.objects.get(
                id=job_h5ad_dataset
            )
            h5ad_dataset_organism = h5ad_dataset_instance.h5ad_dataset_organism.id
            job_annotate_and_plot_data_type = "h5ad"
        except Exception as e:
            job_h5ad_dataset = None
            h5ad_dataset_instance = None

        try:
            tenxfbcm_dataset_instance = TenxfbcmDatasetModel.objects.get(
                id=job_tenxfbcm_dataset
            )
            tenxfbcm_dataset_organism = tenxfbcm_dataset_instance.tenxfbcm_dataset_organism.id
            job_annotate_and_plot_data_type = "10x"
        except Exception as e:
            job_tenxfbcm_dataset = None
            tenxfbcm_dataset_instance = None

        try:
            pred_dataset_instance = PredDatasetModel.objects.get(
                id=job_pred_dataset
            )
        except Exception as e:
            job_pred_dataset = None
            pred_dataset_instance = None

        try:
            predictor_instance = PredictorModel.objects.get(id=job_predictor)
            predictor_organism = predictor_instance.predictor_organism.id
        except Exception as e:
            return Response({"isJobAnnotateAndPlot": False, "error": "Predictor does not exists"}, status=404)

        try:
            if (h5ad_dataset_organism == predictor_organism) or (tenxfbcm_dataset_organism == predictor_organism):
                job_organism = predictor_organism

                pred_dataset_model_serializer = PredDatasetModelSerializer(data={
                    'pred_dataset_name': str(job_name) + "_" + str(job_annotate_and_plot_prediction_filename),
                    'pred_dataset_file_extension': "csv",
                    'pred_dataset_file': ContentFile("\n", name=str(job_annotate_and_plot_prediction_filename) + ".csv"),
                    'pred_dataset_organism': job_organism,
                    'pred_dataset_public_flag': False,
                    'pred_dataset_upload_user': pred_dataset_upload_user
                })

                pred_dataset_model_serializer_instance = None
                if pred_dataset_model_serializer.is_valid():
                    try:
                        pred_dataset_model_serializer_instance = pred_dataset_model_serializer.save()
                    except Exception as e:
                        return Response({"isJobAnnotateAndPlot": False, "error": str(e)}, status=405)
                else:
                    return Response({
                        "isJobAnnotateAndPlot": False,
                        "error": str(pred_dataset_model_serializer.errors)
                    }, status=405)

                job_annotate_and_plot_file_output_model_serializer = JobAnnotateAndPlotFileOutputModelSerializer(
                    data={
                        "job_organism": job_organism,
                        "job_annotate_and_plot_folder": job_annotate_and_plot_folder,
                        "job_annotate_and_plot_output_with_celltype_file": ContentFile("\n", name=str(job_annotate_and_plot_output_with_celltype_filename) + ".h5ad"),
                        "job_annotate_and_plot_prediction_file": ContentFile("\n", name=str(job_annotate_and_plot_prediction_filename) + ".csv"),
                        "job_annotate_and_plot_log_file": ContentFile("\n", name=str(job_annotate_and_plot_log_filename) + ".txt"),
                        "job_annotate_and_plot_stats_csv_file": ContentFile("\n", name=str(job_annotate_and_plot_stats_csv_filename) + ".csv"),
                        "job_annotate_and_plot_stats_pdf_file": ContentFile("\n", name=str(job_annotate_and_plot_stats_pdf_filename) + ".pdf"),
                        "job_annotate_and_plot_marker_genes_file": ContentFile("\n", name=str(job_annotate_and_plot_marker_genes_filename) + ".csv"),
                        "job_annotate_and_plot_annotate_tsne_file": ContentFile("\n", name=str(job_annotate_and_plot_annotate_tsne_filename) + ".pdf"),
                        "job_annotate_and_plot_annotate_umap_file": ContentFile("\n", name=str(job_annotate_and_plot_annotate_umap_filename) + ".pdf"),
                        "job_annotate_and_plot_top3_genes_dotplot_file": ContentFile("\n", name=str(job_annotate_and_plot_top3_genes_dotplot_filename) + ".pdf"),
                        "job_annotate_and_plot_all_markers_file": ContentFile("\n", name=str(job_annotate_and_plot_all_markers_filename) + ".xlsx"),
                        "job_annotate_and_plot_top5_markers_file": ContentFile("\n", name=str(job_annotate_and_plot_top5_markers_filename) + ".xlsx"),
                        "job_annotate_and_plot_top10_markers_file": ContentFile("\n", name=str(job_annotate_and_plot_top10_markers_filename) + ".xlsx"),
                        "job_annotate_and_plot_top25_markers_file": ContentFile("\n", name=str(job_annotate_and_plot_top25_markers_filename) + ".xlsx"),
                        "job_annotate_and_plot_stdout_file": ContentFile("\n", name=str(job_annotate_and_plot_stdout_filename) + ".txt"),
                        "job_annotate_and_plot_stderr_file": ContentFile("\n", name=str(job_annotate_and_plot_stderr_filename) + ".txt"),
                        "job_annotate_and_plot_file_creation_user": job_annotate_and_plot_file_creation_user
                    }
                )

                job_annotate_and_plot_file_output_model_serializer_instance = None
                if job_annotate_and_plot_file_output_model_serializer.is_valid():
                    try:
                        job_annotate_and_plot_file_output_model_serializer_instance = job_annotate_and_plot_file_output_model_serializer.save()
                    except Exception as e:
                        return Response({"isJobAnnotateAndPlot": False, "error": str(e)}, status=405)
                else:
                    return Response({
                        "isJobAnnotateAndPlot": False,
                        "error": str(job_annotate_and_plot_file_output_model_serializer.errors)
                    }, status=405)

                if pred_dataset_model_serializer_instance and job_annotate_and_plot_file_output_model_serializer_instance:

                    job_annotate_and_plot_file_output = job_annotate_and_plot_file_output_model_serializer_instance.id

                    arg_array = [script_instance.script_file.path]

                    if h5ad_dataset_instance:
                        arg_array.append(
                            h5ad_dataset_instance.h5ad_dataset_file.path)
                    else:
                        arg_array.append(None)

                    if tenxfbcm_dataset_instance:
                        arg_array.append(
                            os.path.join(
                                settings.MEDIA_ROOT,
                                tenxfbcm_dataset_instance.tenxfbcm_dataset_folder_path
                            )
                        )
                    else:
                        arg_array.append(None)

                    if pred_dataset_instance:
                        arg_array.append(
                            pred_dataset_instance.pred_dataset_file.path)
                    else:
                        arg_array.append(None)

                    arg_array = arg_array + [
                        predictor_instance.predictor_file.path,
                        job_annotate_and_plot_data_type,
                        os.path.join(
                            settings.MEDIA_ROOT,
                            job_annotate_and_plot_file_output_model_serializer_instance.job_annotate_and_plot_folder_path
                        ),
                        job_annotate_and_plot_file_output_model_serializer_instance.job_annotate_and_plot_prediction_file.path,
                        job_annotate_and_plot_file_output_model_serializer_instance.job_annotate_and_plot_log_file.path,
                        job_annotate_and_plot_file_output_model_serializer_instance.job_annotate_and_plot_stats_csv_file.path,
                        job_annotate_and_plot_file_output_model_serializer_instance.job_annotate_and_plot_stdout_file.path,
                        job_annotate_and_plot_file_output_model_serializer_instance.job_annotate_and_plot_stderr_file.path,
                        pred_dataset_model_serializer_instance.pred_dataset_file.path
                    ]

                    async_result_object = AnnotateAndPlot.apply_async(
                        args=arg_array
                    )

                    job_celery_task_id = async_result_object.id

                    async_result_instance = AsyncResult(job_celery_task_id)

                    job_annotate_and_plot_model_serializer = JobAnnotateAndPlotModelSerializer(
                        data={
                            "job_name": job_name,
                            "job_script": job_script,
                            "job_h5ad_dataset": job_h5ad_dataset,
                            "job_tenxfbcm_dataset": job_tenxfbcm_dataset,
                            "job_pred_dataset": job_pred_dataset,
                            "job_predictor": job_predictor,
                            "job_organism": job_organism,
                            "job_annotate_and_plot_data_type": job_annotate_and_plot_data_type,
                            "job_annotate_and_plot_file_output": job_annotate_and_plot_file_output,
                            "job_celery_task_id": job_celery_task_id,
                            "job_celery_task_status": async_result_instance.status,
                            "job_creation_user": job_creation_user
                        }
                    )

                    if job_annotate_and_plot_model_serializer.is_valid():
                        job_annotate_and_plot_model_serializer_instance = job_annotate_and_plot_model_serializer.save()

                        job_annotate_and_plot_file_output_model_serializer = JobAnnotateAndPlotFileOutputModelSerializer(
                            instance=job_annotate_and_plot_file_output_model_serializer_instance,
                            data={
                                "job_celery_task_id": job_celery_task_id
                            },
                            partial=True
                        )

                        if job_annotate_and_plot_file_output_model_serializer.is_valid():
                            try:
                                job_annotate_and_plot_file_output_model_serializer_instance = job_annotate_and_plot_file_output_model_serializer.save()
                            except Exception as e:
                                try:
                                    if job_annotate_and_plot_file_output_model_serializer_instance:
                                        job_annotate_and_plot_file_output_model_serializer_instance.delete()
                                except Exception as e:
                                    pass
                                try:
                                    if pred_dataset_model_serializer_instance:
                                        pred_dataset_model_serializer_instance.delete()
                                except Exception as e:
                                    pass
                                return Response({"isJobAnnotateAndPlot": False, "error": str(e)}, status=405)
                        else:
                            try:
                                if job_annotate_and_plot_file_output_model_serializer_instance:
                                    job_annotate_and_plot_file_output_model_serializer_instance.delete()
                            except Exception as e:
                                pass
                            try:
                                if pred_dataset_model_serializer_instance:
                                    pred_dataset_model_serializer_instance.delete()
                            except Exception as e:
                                pass
                            return Response({
                                "isJobAnnotateAndPlot": False,
                                "error": str(job_annotate_and_plot_file_output_model_serializer.errors)
                            }, status=405
                            )

                        return Response({"isJobAnnotateAndPlot": True}, status=201)
                    else:
                        try:
                            if job_annotate_and_plot_file_output_model_serializer_instance:
                                job_annotate_and_plot_file_output_model_serializer_instance.delete()
                        except Exception as e:
                            return Response({"isJobAnnotateAndPlot": False, "error": str(e)}, status=405)
                        try:
                            if pred_dataset_model_serializer_instance:
                                pred_dataset_model_serializer_instance.delete()
                        except Exception as e:
                            return Response({"isJobAnnotateAndPlot": False, "error": str(e)}, status=405)
                        return Response({
                            "isJobAnnotateAndPlot": False,
                            "error": str(job_annotate_and_plot_model_serializer.errors)
                        }, status=405
                        )
                else:
                    try:
                        if job_annotate_and_plot_file_output_model_serializer_instance:
                            job_annotate_and_plot_file_output_model_serializer_instance.delete()
                    except Exception as e:
                        pass
                    try:
                        if pred_dataset_model_serializer_instance:
                            pred_dataset_model_serializer_instance.delete()
                    except Exception as e:
                        pass
                    return Response({
                        "isJobAnnotateAndPlot": False,
                        "error": "Arguments and file outputs cannot be created"
                    }, status=405
                    )
            else:
                return Response({
                    "isJobAnnotateAndPlot": False,
                    "error": "Organism of the dataset and predictor do not match"
                }, status=405)
        except Exception as e:
            return Response({"isJobAnnotateAndPlot": False, "error": str(e)}, status=405)

    return Response({"isJobAnnotateAndPlot": False, "error": "Invalid request method"}, status=405)
