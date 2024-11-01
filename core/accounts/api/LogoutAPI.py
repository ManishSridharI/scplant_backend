from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from django.contrib.auth import logout


@api_view(['POST'])
@authentication_classes([])
@permission_classes([])
def Logout(request):
    try:
        logout(request)
        return Response({
            "isLogout": True
        })
    except Exception as e:
        return Response({
            "isLogout": False
        })

