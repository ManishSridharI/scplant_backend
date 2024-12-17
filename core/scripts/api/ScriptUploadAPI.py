from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.ScriptModel import ScriptModel

from ..serializers.ScriptModelSerializer import ScriptModelSerializer


@api_view(['POST'])
def ScriptUpload(request):
    if request.method == 'POST':
        script_file = request.FILES.get('script_file')
        script_name = request.POST.get('script_name')
        script_public_flag = request.POST.get('script_public_flag')
        script_upload_user = request.user.id

        serializer = ScriptModelSerializer(data={
            'script_name': script_name,
            'script_file': script_file,
            'script_public_flag': script_public_flag,
            'script_upload_user': script_upload_user
        })

        # Validate serializer and save if valid
        if serializer.is_valid():
            serializer.save()
            return Response({"isScriptUpload": True}, status=201)
        else:
            return Response({"isScriptUpload": False, "error": str(serializer.errors)}, status=405)

    return Response({"isScriptUpload": False, "error": "Invalid request method"}, status=405)
