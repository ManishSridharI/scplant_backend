import json
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes
from django.contrib.auth import authenticate, login
from ..serializers.CustomUserModelSerializer import CustomUserModelSerializer

@api_view(['POST'])
@authentication_classes([])  # No authentication required for login
@permission_classes([])       # No permission required for login
def Login(request):
    if request.method == 'POST':
        try:
            request_body_dict = json.loads(request.body)
            username = request_body_dict.get('username')  # Fetching username from request
            password = request_body_dict.get('password')  # Fetching password from request

            # Authenticate the user
            user = authenticate(request, username=username, password=password)

            if user is not None:
                login(request, user)  # Log the user in
                response_object = {
                    "isLogin": True,
                    "User": CustomUserModelSerializer(user).data
                }
                return Response(response_object)

        except json.JSONDecodeError:
            return Response({"isLogin": False, "error": "Invalid JSON."}, status=400)
        except Exception as e:
            return Response({"isLogin": False, "error": str(e)}, status=500)

    return Response({"isLogin": False}, status=401)  # Return 401 for unauthorized
