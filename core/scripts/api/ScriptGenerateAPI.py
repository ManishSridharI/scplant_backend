from django.core.files.base import ContentFile, File

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.ScriptModel import ScriptModel

from ..serializers.ScriptModelSerializer import ScriptModelSerializer


@api_view(['POST'])
def ScriptGenerate(request):
    if request.method == 'POST':
        script_name = request.data['script_name']
        script_filename = request.data['script_filename']
        script_file_extension = request.data['script_file_extension']
        script_public_flag = request.data['script_public_flag']
        script_upload_user = request.user.id

        if (script_file_extension == "R") or (script_file_extension == "py"):
            script_model_serializer = ScriptModelSerializer(data={
                'script_name': script_name,
                'script_file_extension': script_file_extension,
                'script_file': ContentFile("\n", name=str(script_filename)+"."+str(script_file_extension)),
                'script_public_flag': script_public_flag,
                'script_upload_user': script_upload_user
            })

            # Validate serializer and save if valid
            if script_model_serializer.is_valid():
                script_model_serializer_instance = script_model_serializer.save()
                return Response({"isScriptGenerate": True, "Script": script_model_serializer_instance.script_file.name}, status=201)
            else:
                return Response({"isScriptGenerate": False, "error": str(script_model_serializer.errors)}, status=405)
        else:
            return Response({"isScriptGenerate": False, "error": "Incorrect script file extension"}, status=405)

    return Response({"isScriptGenerate": False, "error": "Invalid request method"}, status=405)
