from django.db.models import Q

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.ScriptModel import ScriptModel

from ..serializers.ScriptModelSerializer import ScriptModelSerializer


@api_view(['GET'])
def ScriptQuery(request):
    if request.method == 'GET':
        try:
            script_upload_user = request.user.id

            script = ScriptModel.objects.filter(
                script_upload_user=script_upload_user
            )

            response_object = {
                "isScriptQuery": True,
                "Script": ScriptModelSerializer(script, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isScriptQuery": False, "error": str(e)}, status=500)

    return Response({"isScriptQuery": False, "error": "Invalid request method"}, status=405)


@api_view(['GET'])
def ScriptQueryPublic(request):
    if request.method == 'GET':
        try:
            script_upload_user = request.user.id

            script = ScriptModel.objects.filter(
                script_public_flag=True
            )

            response_object = {
                "isScriptQueryPublic": True,
                "Script": ScriptModelSerializer(script, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isScriptQueryPublic": False, "error": str(e)}, status=500)

    return Response({"isScriptQueryPublic": False, "error": "Invalid request method"}, status=405)


@api_view(['GET'])
def ScriptQueryUploadedAndPublic(request):
    if request.method == 'GET':
        try:
            script_upload_user = request.user.id

            script = ScriptModel.objects.filter(
                Q(script_upload_user=script_upload_user) | Q(script_public_flag=True)
            )

            response_object = {
                "isScriptQueryUploadedAndPublic": True,
                "Script": ScriptModelSerializer(script, many=True).data
            }
            return Response(response_object)
        except Exception as e:
            return Response({"isScriptQueryUploadedAndPublic": False, "error": str(e)}, status=500)

    return Response({"isScriptQueryUploadedAndPublic": False, "error": "Invalid request method"}, status=405)
