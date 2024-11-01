import re

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from django.contrib.auth import authenticate, login
from django.contrib.auth.hashers import make_password

from ..models.CustomUserModel import CustomUserModel

from ..serializers.CustomUserModelSerializer import CustomUserModelSerializer


@api_view(['POST'])
@authentication_classes([])
@permission_classes([])
def PasswordChange(request):
    if request.method == 'POST':
        try:
            username = request.POST.get('username')
            password = request.POST.get('password')
            password1 = request.POST.get('password1')
            password2 = request.POST.get('password2')

            user = authenticate(
                request, 
                username = username, 
                password = password
            )

            if user is not None:
                if (password != password1) and (password != password2) and (password1 == password2):
                    if len(password1) > 8:
                        if (re.sub("[a-z]", "", password1) != "") and (re.sub("[A-Z]", "", password1) != "") and (re.sub("[0-9]", "", password1) != ""):
                            user.set_password(password1)

                            try:
                                user.save()
                                response_object = {
                                    "isPasswordChange": True,
                                    "User": CustomUserModelSerializer(user).data
                                }
                                return Response(response_object)
                            except Exception as e:
                                response_object = {
                                    "message": "KKKKK " + str(e),
                                    "isPasswordChange": False
                                }
                                return Response(response_object)

        except Exception as e:
            response_object = {
                "message": "PPPPP " + str(e),
                "isPasswordChange": False
            }
            return Response(response_object)

    response_object = {
        "isPasswordChange": False
    }
    return Response(response_object)

