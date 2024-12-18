from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from celery.result import AsyncResult

from ..serializers.JobAnnotateAndPlotModelSerializer import JobAnnotateAndPlotModelSerializer
from ..serializers.JobAnnotateAndPlotArgumentModelSerializer import JobAnnotateAndPlotArgumentModelSerializer
from ..serializers.JobAnnotateAndPlotFileOutputModelSerializer import JobAnnotateAndPlotFileOutputModelSerializer

from ..tasks.AnnotateAndPlotTask import AnnotateAndPlot

from scripts.models.ScriptModel import ScriptModel
from datasets.models.DatasetModel import DatasetModel
from predictors.models.PredictorModel import PredictorModel


@api_view(['POST'])
def JobAnnotateAndPlot(request):
    if request.method == 'POST':
        job_name = request.data['job_name']
        job_script = request.data['job_script']
        job_dataset = request.data['job_dataset']
        job_predictor = request.data['job_predictor']
        job_annotate_and_plot_gene_number = request.data['job_annotate_and_plot_gene_number']
        job_annotate_and_plot_log_filename = request.data['job_annotate_and_plot_log_filename']
        job_annotate_and_plot_stdout_filename = request.data['job_annotate_and_plot_stdout_filename']
        job_annotate_and_plot_stderr_filename = request.data['job_annotate_and_plot_stderr_filename']

        job_annotate_and_plot_top25_markers_filename = "top25_markers"
        job_annotate_and_plot_marker_genes_filename = "marker_genes"
        job_annotate_and_plot_output_with_celltype_filename = "output_with_celltype"
        job_annotate_and_plot_prediction_filename = "prediction"
        job_annotate_and_plot_annotate_tsne_filename = "annotate_tsne"
        job_annotate_and_plot_annotate_umap_filename = "annotate_umap"
        job_annotate_and_plot_top3_genes_dotplot_filename = "top3_genes_dotplot"

        job_creation_user = request.user.id
        job_annotate_and_plot_file_creation_user = request.user.id
        job_annotate_and_plot_argument_creation_user = request.user.id

        try:
            script_instance = ScriptModel.objects.get(id=job_script)
        except Exception as e:
            return Response({"isJobAnnotateAndPlot": False, "error": "Script does not exists"}, status=404)

        try:
            dataset_instance = DatasetModel.objects.get(id=job_dataset)
        except Exception as e:
            return Response({"isJobAnnotateAndPlot": False, "error": "Dataset does not exists"}, status=404)

        try:
            predictor_instance = PredictorModel.objects.get(id=job_predictor)
        except Exception as e:
            return Response({"isJobAnnotateAndPlot": False, "error": "Predictor does not exists"}, status=404)

        try:
            job_annotate_and_plot_argument_model_serializer = JobAnnotateAndPlotArgumentModelSerializer(
                data={
                    "job_annotate_and_plot_gene_number": job_annotate_and_plot_gene_number,
                    "job_annotate_and_plot_argument_creation_user": job_annotate_and_plot_argument_creation_user
                }
            )

            job_annotate_and_plot_argument_model_serializer_instance = None
            if job_annotate_and_plot_argument_model_serializer.is_valid():
                try:
                    job_annotate_and_plot_argument_model_serializer_instance = job_annotate_and_plot_argument_model_serializer.save()
                except Exception as e:
                    return Response({"isJobAnnotateAndPlot": False, "error": str(e)}, status=405)
            else:
                return Response({"isJobAnnotateAndPlot": False, "error": str(job_annotate_and_plot_argument_model_serializer.errors)}, status=405)

            job_annotate_and_plot_file_output_model_serializer = JobAnnotateAndPlotFileOutputModelSerializer(
                data={
                    "job_annotate_and_plot_top25_markers_file": ContentFile("\n", name=str(job_annotate_and_plot_top25_markers_filename)+".txt"),
                    "job_annotate_and_plot_marker_genes_file": ContentFile("\n", name=str(job_annotate_and_plot_marker_genes_filename)+".csv"),
                    "job_annotate_and_plot_output_with_celltype_file": ContentFile("\n", name=str(job_annotate_and_plot_output_with_celltype_filename)+".h5ad"),
                    "job_annotate_and_plot_log_file": ContentFile("\n", name=str(job_annotate_and_plot_log_filename)+".txt"),
                    "job_annotate_and_plot_prediction_file": ContentFile("\n", name=str(job_annotate_and_plot_prediction_filename)+".csv"),
                    "job_annotate_and_plot_annotate_tsne_file": ContentFile("\n", name=str(job_annotate_and_plot_annotate_tsne_filename)+".pdf"),
                    "job_annotate_and_plot_annotate_umap_file": ContentFile("\n", name=str(job_annotate_and_plot_annotate_umap_filename)+".pdf"),
                    "job_annotate_and_plot_top3_genes_dotplot_file": ContentFile("\n", name=str(job_annotate_and_plot_top3_genes_dotplot_filename)+".pdf"),
                    "job_annotate_and_plot_stdout_file": ContentFile("\n", name=str(job_annotate_and_plot_stdout_filename)+".txt"),
                    "job_annotate_and_plot_stderr_file": ContentFile("\n", name=str(job_annotate_and_plot_stderr_filename)+".txt"),
                    "job_annotate_and_plot_file_creation_user": job_annotate_and_plot_file_creation_user
                }
            )

            job_annotate_and_plot_file_output_model_serializer_instance = None
            if job_annotate_and_plot_file_output_model_serializer.is_valid():
                try:
                    job_annotate_and_plot_file_output_model_serializer_instance = job_annotate_and_plot_file_output_model_serializer.save()
                except Exception as e:
                    try:
                        if job_annotate_and_plot_argument_model_serializer_instance:
                            job_annotate_and_plot_argument_model_serializer_instance.delete()
                    except Exception as e:
                        return Response({"isJobAnnotateAndPlot": False, "error": str(e)}, status=405)
                    return Response({"isJobAnnotateAndPlot": False, "error": str(e)}, status=405)
            else:
                try:
                    if job_annotate_and_plot_argument_model_serializer_instance:
                        job_annotate_and_plot_argument_model_serializer_instance.delete()
                except Exception as e:
                    return Response({"isJobAnnotateAndPlot": False, "error": str(e)}, status=405)
                return Response({"isJobAnnotateAndPlot": False, "error": str(job_annotate_and_plot_file_output_model_serializer.errors)}, status=405)

            if job_annotate_and_plot_argument_model_serializer_instance and job_annotate_and_plot_file_output_model_serializer_instance:

                job_annotate_and_plot_argument = job_annotate_and_plot_argument_model_serializer_instance.id
                job_annotate_and_plot_file_output = job_annotate_and_plot_file_output_model_serializer_instance.id

                async_result_object = AnnotateAndPlot.apply_async(
                    args=[
                        script_instance.script_file.path,
                        dataset_instance.dataset_file.path,
                        predictor_instance.predictor_file.path,
                        job_annotate_and_plot_argument_model_serializer_instance.job_annotate_and_plot_gene_number,
                        job_annotate_and_plot_file_output_model_serializer_instance.job_annotate_and_plot_top25_markers_file.path,
                        job_annotate_and_plot_file_output_model_serializer_instance.job_annotate_and_plot_marker_genes_file.path,
                        job_annotate_and_plot_file_output_model_serializer_instance.job_annotate_and_plot_output_with_celltype_file.path,
                        job_annotate_and_plot_file_output_model_serializer_instance.job_annotate_and_plot_log_file.path,
                        job_annotate_and_plot_file_output_model_serializer_instance.job_annotate_and_plot_prediction_file.path,
                        job_annotate_and_plot_file_output_model_serializer_instance.job_annotate_and_plot_annotate_tsne_file.path,
                        job_annotate_and_plot_file_output_model_serializer_instance.job_annotate_and_plot_annotate_umap_file.path,
                        job_annotate_and_plot_file_output_model_serializer_instance.job_annotate_and_plot_top3_genes_dotplot_file.path,
                        job_annotate_and_plot_file_output_model_serializer_instance.job_annotate_and_plot_stdout_file.path,
                        job_annotate_and_plot_file_output_model_serializer_instance.job_annotate_and_plot_stderr_file.path
                    ]
                )

                job_celery_task_id = async_result_object.id

                async_result_instance = AsyncResult(job_celery_task_id)

                job_annotate_and_plot_model_serializer = JobAnnotateAndPlotModelSerializer(
                    data={
                        "job_name": job_name,
                        "job_script": job_script,
                        "job_dataset": job_dataset,
                        "job_predictor": job_predictor,
                        "job_annotate_and_plot_argument": job_annotate_and_plot_argument,
                        "job_annotate_and_plot_file_output": job_annotate_and_plot_file_output,
                        "job_celery_task_id": job_celery_task_id,
                        "job_celery_task_status": async_result_instance.status,
                        "job_creation_user": job_creation_user
                    }
                )

                if job_annotate_and_plot_model_serializer.is_valid():
                    job_annotate_and_plot_model_serializer_instance = job_annotate_and_plot_model_serializer.save()

                    job_annotate_and_plot_argument_model_serializer = JobAnnotateAndPlotArgumentModelSerializer(
                        instance=job_annotate_and_plot_argument_model_serializer_instance,
                        data={
                            "job_celery_task_id": job_celery_task_id
                        },
                        partial=True
                    )

                    if job_annotate_and_plot_argument_model_serializer.is_valid():
                        try:
                            job_annotate_and_plot_argument_model_serializer_instance = job_annotate_and_plot_argument_model_serializer.save()
                        except Exception as e:
                            return Response({"isJobAnnotateAndPlot": False, "error": str(e)}, status=405)
                    else:
                        return Response({"isJobAnnotateAndPlot": False, "error": str(job_annotate_and_plot_argument_model_serializer.errors)}, status=405)

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
                            return Response({"isJobAnnotateAndPlot": False, "error": str(e)}, status=405)
                    else:
                        return Response({"isJobAnnotateAndPlot": False, "error": str(job_annotate_and_plot_file_output_model_serializer.errors)}, status=405)

                    return Response({"isJobAnnotateAndPlot": True}, status=201)
                else:
                    try:
                        if job_annotate_and_plot_argument_model_serializer_instance:
                            job_annotate_and_plot_argument_model_serializer_instance.delete()
                    except Exception as e:
                        return Response({"isJobAnnotateAndPlot": False, "error": str(e)}, status=405)
                    try:
                        if job_annotate_and_plot_file_output_model_serializer_instance:
                            job_annotate_and_plot_file_output_model_serializer_instance.delete()
                    except Exception as e:
                        return Response({"isJobAnnotateAndPlot": False, "error": str(e)}, status=405)
                    return Response({"isJobAnnotateAndPlot": False, "error": str(job_annotate_and_plot_model_serializer.errors)}, status=405)
            else:
                try:
                    if job_annotate_and_plot_argument_model_serializer_instance:
                        job_annotate_and_plot_argument_model_serializer_instance.delete()
                except Exception as e:
                    pass
                try:
                    if job_annotate_and_plot_file_output_model_serializer_instance:
                        job_annotate_and_plot_file_output_model_serializer_instance.delete()
                except Exception as e:
                    pass
                return Response({"isJobAnnotateAndPlot": False, "error": str(job_annotate_and_plot_argument_model_serializer.errors)}, status=405)

        except Exception as e:
            return Response({"isJobAnnotateAndPlot": False, "error": str(e)}, status=405)

    return Response({"isJobAnnotateAndPlot": False, "error": "Invalid request method"}, status=405)
