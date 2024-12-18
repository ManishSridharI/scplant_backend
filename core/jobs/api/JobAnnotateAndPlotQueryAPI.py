from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.JobAnnotateAndPlotModel import JobAnnotateAndPlotModel

from ..serializers.JobAnnotateAndPlotModelSerializer import JobAnnotateAndPlotModelSerializer


@api_view(['GET'])
def JobAnnotateAndPlotQuery(request):
    if request.method == 'GET':
        try:
            job_creation_user = request.user.id

            job_annotate_and_plot = JobAnnotateAndPlotModel.objects.filter(
                job_creation_user=job_creation_user
            )

            response_object = {
                "isJobAnnotateAndPlotQuery": True,
                "JobAnnotateAndPlot": JobAnnotateAndPlotModelSerializer(job_annotate_and_plot, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isJobAnnotateAndPlotQuery": False, "error": str(e)}, status=500)

    return Response({"isJobAnnotateAndPlotQuery": False, "error": "Invalid request method"}, status=405)
