from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.JobTreatmentVsControlModel import JobTreatmentVsControlModel

from ..serializers.JobTreatmentVsControlModelSerializer import JobTreatmentVsControlModelSerializer


@api_view(['GET'])
def JobTreatmentVsControlQuery(request):
    if request.method == 'GET':
        try:
            job_creation_user = request.user.id

            job_treatment_vs_control = JobTreatmentVsControlModel.objects.filter(
                job_creation_user=job_creation_user
            )

            response_object = {
                "isJobTreatmentVsControlQuery": True,
                "JobTreatmentVsControl": JobTreatmentVsControlModelSerializer(job_treatment_vs_control, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isJobTreatmentVsControlQuery": False, "error": str(e)}, status=500)

    return Response({"isJobTreatmentVsControlQuery": False, "error": "Invalid request method"}, status=405)
