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
            username = request.data['username']
            email = request.data['email']
            password1 = request.data['password1']
            password2 = request.data['password2']
            first_name = request.data['first_name']
            last_name = request.data['last_name']
            organization = request.data['organization']

            # Validate inputs
            if not all([username, email, password1, password2]):
                return Response({"isRegister": False, "error": "Missing required fields"}, status=400)

            if password1 != password2:
                return Response({"isRegister": False, "error": "Passwords do not match"}, status=400)

            if len(password1) <= 8:
                return Response({"isRegister": False, "error": "Password must be longer than 8 characters"}, status=400)

            if not re.search("[a-z]", password1) or not re.search("[A-Z]", password1) or not re.search("[0-9]", password1):
                return Response({"isRegister": False, "error": "Password must contain uppercase, lowercase, and digits"}, status=400)

            # Create user serializer
            serializer = CustomUserModelSerializer(data={
                'username': username,
                'email': email,
                'password': make_password(password1),
                'first_name': first_name,
                'last_name': last_name,
                'organization': organization
            })

            # Validate serializer and save if valid
            if serializer.is_valid():
                serializer.save()
                return Response({"isRegister": True}, status=201)
            else:
                # If serializer is not valid, return errors
                return Response({"isRegister": False, "error": serializer.errors}, status=400)

        except Exception as e:
            return Response({"isRegister": False, "error": str(e)}, status=500)

    return Response({"isRegister": False, "error": "Invalid request method"}, status=405)
