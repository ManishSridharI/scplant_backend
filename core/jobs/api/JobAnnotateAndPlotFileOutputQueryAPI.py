from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.JobAnnotateAndPlotFileOutputModel import JobAnnotateAndPlotFileOutputModel

from ..serializers.JobAnnotateAndPlotFileOutputModelSerializer import JobAnnotateAndPlotFileOutputModelSerializer


@api_view(['GET'])
def JobAnnotateAndPlotFileOutputQuery(request):
    if request.method == 'GET':
        try:
            job_annotate_and_plot_file_creation_user = request.user.id

            job_annotate_and_plot_file_output_instance = JobAnnotateAndPlotFileOutputModel.objects.filter(
                job_annotate_and_plot_file_creation_user=job_annotate_and_plot_file_creation_user
            )

            response_object = {
                "isJobAnnotateAndPlotFileOutputQuery": True,
                "JobAnnotateAndPlotFileOutput": JobAnnotateAndPlotFileOutputModelSerializer(job_annotate_and_plot_file_output_instance, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isJobAnnotateAndPlotFileOutputQuery": False, "error": str(e)}, status=500)

    return Response({"isJobAnnotateAndPlotFileOutputQuery": False, "error": "Invalid request method"}, status=405)


@api_view(['POST'])
def JobAnnotateAndPlotFileOutputQueryByID(request):
    if request.method == 'POST':
        try:
            job_annotate_and_plot_file_output_id = request.data['job_annotate_and_plot_file_output_id']
            job_annotate_and_plot_file_creation_user = request.user.id

            job_annotate_and_plot_file_output_instance = JobAnnotateAndPlotFileOutputModel.objects.get(
                pk=job_annotate_and_plot_file_output_id,
                job_annotate_and_plot_file_creation_user=job_annotate_and_plot_file_creation_user
            )

            response_object = {
                "isJobAnnotateAndPlotFileOutputQueryByID": True,
                "JobAnnotateAndPlotFileOutput": JobAnnotateAndPlotFileOutputModelSerializer(job_annotate_and_plot_file_output_instance).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isJobAnnotateAndPlotFileOutputQueryByID": False, "error": str(e)}, status=500)

    return Response({"isJobAnnotateAndPlotFileOutputQueryByID": False, "error": "Invalid request method"}, status=405)
