from rest_framework.decorators import api_view
from rest_framework.response import Response

from ..models.TestModel import TestModel

from ..serializers.TestModelSerializer import TestModelSerializer


@api_view(['GET'])
def Test(request):
    if True:
        tests = TestModel.objects.all()
        return Response(TestModelSerializer(tests, many=True).data)
    else:
        return Response({"Test": "Test"})

