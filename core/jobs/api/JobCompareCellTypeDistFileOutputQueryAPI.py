from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.JobCompareCellTypeDistFileOutputModel import JobCompareCellTypeDistFileOutputModel

from ..serializers.JobCompareCellTypeDistFileOutputModelSerializer import JobCompareCellTypeDistFileOutputModelSerializer


@api_view(['GET'])
def JobCompareCellTypeDistFileOutputQuery(request):
    if request.method == 'GET':
        try:
            job_compare_cell_type_dist_file_creation_user = request.user.id

            job_compare_cell_type_dist_file_output = JobCompareCellTypeDistFileOutputModel.objects.filter(
                job_compare_cell_type_dist_file_creation_user=job_compare_cell_type_dist_file_creation_user
            )

            response_object = {
                "isJobCompareCellTypeDistFileOutputQuery": True,
                "JobCompareCellTypeDistFileOutput": JobCompareCellTypeDistFileOutputModelSerializer(job_compare_cell_type_dist_file_output, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isJobCompareCellTypeDistFileOutputQuery": False, "error": str(e)}, status=500)

    return Response({"isJobCompareCellTypeDistFileOutputQuery": False, "error": "Invalid request method"}, status=405)


@api_view(['POST'])
def JobCompareCellTypeDistFileOutputQueryByID(request):
    if request.method == 'POST':
        try:
            job_compare_cell_type_dist_file_output_id = request.data['job_compare_cell_type_dist_file_output_id']
            job_compare_cell_type_dist_file_creation_user = request.user.id

            job_compare_cell_type_dist_file_output = JobCompareCellTypeDistFileOutputModel.objects.get(
                pk=job_compare_cell_type_dist_file_output_id,
                job_compare_cell_type_dist_file_creation_user=job_compare_cell_type_dist_file_creation_user
            )

            response_object = {
                "isJobCompareCellTypeDistFileOutputQueryByID": True,
                "JobCompareCellTypeDistFileOutput": JobCompareCellTypeDistFileOutputModelSerializer(job_compare_cell_type_dist_file_output).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isJobCompareCellTypeDistFileOutputQueryByID": False, "error": str(e)}, status=500)

    return Response({"isJobCompareCellTypeDistFileOutputQueryByID": False, "error": "Invalid request method"}, status=405)
