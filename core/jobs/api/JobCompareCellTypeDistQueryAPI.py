from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.JobCompareCellTypeDistModel import JobCompareCellTypeDistModel

from ..serializers.JobCompareCellTypeDistModelSerializer import JobCompareCellTypeDistModelSerializer


@api_view(['GET'])
def JobCompareCellTypeDistQuery(request):
    if request.method == 'GET':
        try:
            job_creation_user = request.user.id

            job_compare_cell_type_dist = JobCompareCellTypeDistModel.objects.filter(
                job_creation_user=job_creation_user
            )

            response_object = {
                "isJobCompareCellTypeDistQuery": True,
                "JobCompareCellTypeDist": JobCompareCellTypeDistModelSerializer(job_compare_cell_type_dist, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isJobCompareCellTypeDistQuery": False, "error": str(e)}, status=500)

    return Response({"isJobCompareCellTypeDistQuery": False, "error": "Invalid request method"}, status=405)
