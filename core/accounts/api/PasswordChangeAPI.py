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
            username = request.data['username']
            password = request.data['password']
            password1 = request.data['password1']
            password2 = request.data['password2']

            # For form only
            # username = request.POST.get('username')
            # password = request.POST.get('password')
            # password1 = request.POST.get('password1')
            # password2 = request.POST.get('password2')

            if (username is not None) and (password is not None) and (password1 is not None) and (password2 is not None):
                if (username != "") and (password != "") and (password1 != "") and (password2 != ""):
                    if (password != password1) and (password != password2) and (password1 == password2):
                        if len(password1) > 8:
                            if (re.sub("[a-z]", "", password1) != "") and (re.sub("[A-Z]", "", password1) != "") and (re.sub("[0-9]", "", password1) != ""):
                                user = authenticate(
                                    request,
                                    username=username,
                                    password=password
                                )

                                if user is not None:
                                    serializer = CustomUserModelSerializer(
                                        instance=user,
                                        data={
                                            'username': user.username,
                                            'email': user.email,
                                            'password': make_password(password1),
                                            'first_name': user.first_name,
                                            'last_name': user.last_name,
                                            'organization': user.organization
                                        }
                                    )

                                    if serializer.is_valid():
                                        try:
                                            serializer.save()
                                            response_object = {
                                                "isPasswordChange": True,
                                                "User": CustomUserModelSerializer(user).data
                                            }
                                            return Response(response_object)
                                        except Exception as e:
                                            response_object = {
                                                "isPasswordChange": False
                                            }
                                            return Response(response_object)
                                    else:
                                        response_object = {
                                            "isPasswordChange": False,
                                            "error": serializer.errors
                                        }
                                        return Response(response_object, status=400)

        except Exception as e:
            response_object = {
                "isPasswordChange": False
            }
            return Response(response_object)

    response_object = {
        "isPasswordChange": False
    }
    return Response(response_object)
