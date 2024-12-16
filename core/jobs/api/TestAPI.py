from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.decorators import authentication_classes, permission_classes

from ..tasks.TestTask import Add, Subtract, Multiply, Divide, WriteDate


@api_view(['GET'])
@authentication_classes([])
@permission_classes([])
def Test(request):
    # async_result_object = Add.delay(4, 5)
    # async_result_object = Subtract.delay(5, 3)
    # async_result_object = Multiply.delay(2, 3)
    # async_result_object = Divide.delay(6, 3)
    async_result_object = WriteDate.delay()
    return Response({"isTest": True}, status=201)
