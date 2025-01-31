from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.JobAnnotateAndPlotModel import JobAnnotateAndPlotModel

from ..serializers.JobAnnotateAndPlotModelSerializer import JobAnnotateAndPlotModelSerializer
from ..serializers.JobAnnotateAndPlotFileOutputModelSerializer import JobAnnotateAndPlotFileOutputModelSerializer

from h5addatasets.serializers.H5adDatasetModelSerializer import H5adDatasetModelSerializer
from tenxfeaturebcmatrixdatasets.serializers.TenxfbcmDatasetModelSerializer import TenxfbcmDatasetModelSerializer
from preddatasets.serializers.PredDatasetModelSerializer import PredDatasetModelSerializer


@api_view(['GET'])
def JobAnnotateAndPlotQuery(request):
    if request.method == 'GET':
        try:
            job_creation_user = request.user.id

            job_annotate_and_plot_instance = JobAnnotateAndPlotModel.objects.filter(
                job_creation_user=job_creation_user
            )

            response_object = {
                "isJobAnnotateAndPlotQuery": True,
                "JobAnnotateAndPlot": JobAnnotateAndPlotModelSerializer(job_annotate_and_plot_instance, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isJobAnnotateAndPlotQuery": False, "error": str(e)}, status=500)

    return Response({"isJobAnnotateAndPlotQuery": False, "error": "Invalid request method"}, status=405)


@api_view(['POST'])
def JobAnnotateAndPlotQueryByID(request):
    if request.method == 'POST':
        try:
            job_annotate_and_plot_id = request.data['job_annotate_and_plot_id']

            job_annotate_and_plot_instance = JobAnnotateAndPlotModel.objects.get(
                id=job_annotate_and_plot_id
            )

            job_annotate_and_plot_model_serializer = JobAnnotateAndPlotModelSerializer(
                job_annotate_and_plot_instance)

            job_annotate_and_plot_file_output_instance = job_annotate_and_plot_instance.job_annotate_and_plot_file_output
            job_annotate_and_plot_file_output_model_serializer = JobAnnotateAndPlotFileOutputModelSerializer(
                job_annotate_and_plot_file_output_instance
            )

            response_object = {
                "isJobAnnotateAndPlotQuery": True,
                "JobAnnotateAndPlot": job_annotate_and_plot_model_serializer.data,
                "JobAnnotateAndPlotFileOutput": job_annotate_and_plot_file_output_model_serializer.data
            }

            if job_annotate_and_plot_instance.job_h5ad_dataset:
                h5ad_dataset_instance = job_annotate_and_plot_instance.job_h5ad_dataset
                h5ad_dataset_model_serializer = H5adDatasetModelSerializer(
                    h5ad_dataset_instance
                )
                response_object["H5adDataset"] = h5ad_dataset_model_serializer.data

            if job_annotate_and_plot_instance.job_tenxfbcm_dataset:
                tenxfbcm_dataset_instance = job_annotate_and_plot_instance.job_tenxfbcm_dataset
                tenxfbcm_dataset_model_serializer = TenxfbcmDatasetModelSerializer(
                    tenxfbcm_dataset_instance
                )
                response_object["TenxfbcmDataset"] = tenxfbcm_dataset_model_serializer.data

            if job_annotate_and_plot_instance.job_pred_dataset:
                pred_dataset_instance = job_annotate_and_plot_instance.job_pred_dataset
                pred_dataset_model_serializer = PredDatasetModelSerializer(
                    pred_dataset_instance
                )
                response_object["PredDataset"] = pred_dataset_model_serializer.data

            return Response(response_object)
        except Exception as e:
            return Response({"isJobAnnotateAndPlotQuery": False, "error": str(e)}, status=500)

    return Response({"isJobAnnotateAndPlotQuery": False, "error": "Invalid request method"}, status=405)
