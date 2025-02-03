from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.JobConvertRdsToH5adModel import JobConvertRdsToH5adModel

from ..serializers.JobConvertRdsToH5adModelSerializer import JobConvertRdsToH5adModelSerializer


@api_view(['GET'])
def JobConvertRdsToH5adQuery(request):
    if request.method == 'GET':
        try:
            job_creation_user = request.user.id

            job_convert_rds_to_h5ad_instance = JobConvertRdsToH5adModel.objects.filter(
                job_creation_user=job_creation_user
            )

            response_object = {
                "isJobConvertRdsToH5adQuery": True,
                "JobConvertRdsToH5ad": JobConvertRdsToH5adModelSerializer(job_convert_rds_to_h5ad_instance, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isJobConvertRdsToH5adQuery": False, "error": str(e)}, status=500)

    return Response({"isJobConvertRdsToH5adQuery": False, "error": "Invalid request method"}, status=405)
