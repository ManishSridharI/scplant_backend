import json

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..models.ScriptModel import ScriptModel

from ..serializers.ScriptModelSerializer import ScriptModelSerializer


@api_view(['POST'])
def ScriptDelete(request):
    if request.method == 'POST':
        try:
            script_id = request.data['script_id']

            script_instance = ScriptModel.objects.filter(id=script_id)
            if script_instance:
                script_instance.delete()
                response_object = {"isScriptDelete": True}
                return Response(response_object)
            else:
                response_object = {
                    "isScriptDelete": False,
                    "error": "Data does not exists"
                }
                return Response(response_object, status=404)
        except Exception as e:
            return Response({"isScriptDelete": False, "error": str(e)}, status=500)

    return Response({"isScriptDelete": False, "error": "Invalid request method"}, status=405)
