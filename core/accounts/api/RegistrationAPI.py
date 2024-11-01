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

            if (username is not None) and (email is not None) and (password1 is not None) and (password2 is not None):
                if (username != "") and (email != "") and (password1 != "") and (password2 != ""):
                    if password1 == password2:
                        if len(password1) > 8:
                            if (re.sub("[a-z]", "", password1) != "") and (re.sub("[A-Z]", "", password1) != "") and (re.sub("[0-9]", "", password1) != ""):
                                serializer = CustomUserModelSerializer(
                                    data={
                                        'username': username,
                                        'email': email,
                                        'password': make_password(password1),
                                        'first_name': first_name,
                                        'last_name': last_name,
                                        'organization': organization
                                    }
                                )

                                if serializer.is_valid():
                                    try:
                                        serializer.save()
                                        response_object = {
                                            "isRegister": True
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

