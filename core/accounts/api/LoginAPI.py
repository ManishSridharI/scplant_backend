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
            username = request.data['username']
            password = request.data['password']

            # Authenticate the user
            user = authenticate(request, username=username, password=password)

            if user is not None:
                login(request, user)  # Log the user in
                response_object = {
                    "isLogin": True,
                    "User": CustomUserModelSerializer(user).data
                }
                return Response(response_object)

        except Exception as e:
            return Response({"isLogin": False, "error": str(e)}, status=500)

    # Return 401 for unauthorized
    return Response({"isLogin": False}, status=401)
