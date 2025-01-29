from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.JobConvertRdsToH5adFileOutputModel import JobConvertRdsToH5adFileOutputModel

from ..serializers.JobConvertRdsToH5adFileOutputModelSerializer import JobConvertRdsToH5adFileOutputModelSerializer


@api_view(['GET'])
def JobConvertRdsToH5adFileOutputQuery(request):
    if request.method == 'GET':
        try:
            job_convert_rds_to_h5ad_file_creation_user = request.user.id

            job_convert_rds_to_h5ad_file_output_instance = JobConvertRdsToH5adFileOutputModel.objects.filter(
                job_convert_rds_to_h5ad_file_creation_user=job_convert_rds_to_h5ad_file_creation_user
            )

            response_object = {
                "isJobConvertRdsToH5adFileOutputQuery": True,
                "JobConvertRdsToH5adFileOutput": JobConvertRdsToH5adFileOutputModelSerializer(job_convert_rds_to_h5ad_file_output_instance, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isJobConvertRdsToH5adFileOutputQuery": False, "error": str(e)}, status=500)

    return Response({"isJobConvertRdsToH5adFileOutputQuery": False, "error": "Invalid request method"}, status=405)


@api_view(['POST'])
def JobConvertRdsToH5adFileOutputQueryByID(request):
    if request.method == 'POST':
        try:
            job_convert_rds_to_h5ad_file_output_id = request.data['job_convert_rds_to_h5ad_file_output_id']
            job_convert_rds_to_h5ad_file_creation_user = request.user.id

            job_convert_rds_to_h5ad_file_output_instance = JobConvertRdsToH5adFileOutputModel.objects.get(
                pk=job_convert_rds_to_h5ad_file_output_id,
                job_convert_rds_to_h5ad_file_creation_user=job_convert_rds_to_h5ad_file_creation_user
            )

            response_object = {
                "isJobConvertRdsToH5adFileOutputQueryByID": True,
                "JobConvertRdsToH5adFileOutput": JobConvertRdsToH5adFileOutputModelSerializer(job_convert_rds_to_h5ad_file_output_instance).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isJobConvertRdsToH5adFileOutputQueryByID": False, "error": str(e)}, status=500)

    return Response({"isJobConvertRdsToH5adFileOutputQueryByID": False, "error": "Invalid request method"}, status=405)
