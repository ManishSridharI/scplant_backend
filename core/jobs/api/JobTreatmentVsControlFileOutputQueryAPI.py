from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.JobTreatmentVsControlFileOutputModel import JobTreatmentVsControlFileOutputModel

from ..serializers.JobTreatmentVsControlFileOutputModelSerializer import JobTreatmentVsControlFileOutputModelSerializer


@api_view(['GET'])
def JobTreatmentVsControlFileOutputQuery(request):
    if request.method == 'GET':
        try:
            job_treatment_vs_control_file_creation_user = request.user.id

            job_treatment_vs_control_file_output_instance = JobTreatmentVsControlFileOutputModel.objects.filter(
                job_treatment_vs_control_file_creation_user=job_treatment_vs_control_file_creation_user
            )

            response_object = {
                "isJobTreatmentVsControlFileOutputQuery": True,
                "JobTreatmentVsControlFileOutput": JobTreatmentVsControlFileOutputModelSerializer(job_treatment_vs_control_file_output_instance, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isJobTreatmentVsControlFileOutputQuery": False, "error": str(e)}, status=500)

    return Response({"isJobTreatmentVsControlFileOutputQuery": False, "error": "Invalid request method"}, status=405)


@api_view(['POST'])
def JobTreatmentVsControlFileOutputQueryByID(request):
    if request.method == 'POST':
        try:
            job_treatment_vs_control_file_output_id = request.data['job_treatment_vs_control_file_output_id']
            job_treatment_vs_control_file_creation_user = request.user.id

            job_treatment_vs_control_file_output_instance = JobTreatmentVsControlFileOutputModel.objects.get(
                pk=job_treatment_vs_control_file_output_id,
                job_treatment_vs_control_file_creation_user=job_treatment_vs_control_file_creation_user
            )

            response_object = {
                "isJobTreatmentVsControlFileOutputQueryByID": True,
                "JobTreatmentVsControlFileOutput": JobTreatmentVsControlFileOutputModelSerializer(job_treatment_vs_control_file_output_instance).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isJobTreatmentVsControlFileOutputQueryByID": False, "error": str(e)}, status=500)

    return Response({"isJobTreatmentVsControlFileOutputQueryByID": False, "error": "Invalid request method"}, status=405)
