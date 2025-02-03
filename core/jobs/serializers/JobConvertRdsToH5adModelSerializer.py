from rest_framework import serializers

from ..models.JobConvertRdsToH5adModel import JobConvertRdsToH5adModel


class JobConvertRdsToH5adModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobConvertRdsToH5adModel
        fields = '__all__'
