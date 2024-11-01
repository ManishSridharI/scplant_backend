import re

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from django.contrib.auth.hashers import make_password

from ..models.CustomUserModel import CustomUserModel

from ..serializers.CustomUserModelSerializer import CustomUserModelSerializer


@api_view(['POST'])
@authentication_classes([])
@permission_classes([])
def Registration(request):
    if request.method == 'POST':
        try:
            username = request.POST.get('username')
            email = request.POST.get('email')
            password1 = request.POST.get('password1')
            password2 = request.POST.get('password2')
            first_name = request.POST.get('first_name')
            last_name = request.POST.get('last_name')
            organization = request.POST.get('organization')

            if (username != "") and (email != ""):
                if ("@" in email) and ("." in email):
                    if password1 == password2:
                        if len(password1) > 8:
                            if (re.sub("[a-z]", "", password1) != "") and (re.sub("[A-Z]", "", password1) != "") and (re.sub("[0-9]", "", password1) != ""):
                                user = CustomUserModel(
                                    username=username,
                                    email=email,
                                    first_name=first_name,
                                    last_name=last_name,
                                    organization=organization,
                                )
                                user.set_password(password1)

                                try:
                                    user.save()
                                    response_object = {
                                        "isRegister": True,
                                        "User": CustomUserModelSerializer(user).data
                                    }
                                    return Response(response_object)
                                except Exception as e:
                                    response_object = {
                                        "isRegister": False
                                    }
                                    return Response(response_object)

        except Exception as e:
            response_object = {
                "isRegister": False
            }
            return Response(response_object)

    response_object = {
        "isRegister": False
    }
    return Response(response_object)

