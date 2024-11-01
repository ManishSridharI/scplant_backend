from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from django.contrib.auth import authenticate, login

from ..serializers.CustomUserModelSerializer import CustomUserModelSerializer


@api_view(['POST'])
@authentication_classes([])
@permission_classes([])
def Login(request):
    if request.method == 'POST':
        try:
            username = request.POST.get('username')
            password = request.POST.get('password')

            user = authenticate(
                request, 
                username = username, 
                password = password
            )

            if user is not None:
                login(request, user)
                response_object = {
                    "isLogin": True,
                    "User": CustomUserModelSerializer(user).data
                }
                return Response(response_object)

        except Exception as e:
            response_object = {
                "isLogin": False
            }
            return Response(response_object)

    response_object = {
        "isLogin": False
    }
    return Response(response_object)

