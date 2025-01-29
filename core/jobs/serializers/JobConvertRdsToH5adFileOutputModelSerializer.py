from rest_framework import serializers

from ..models.JobConvertRdsToH5adFileOutputModel import JobConvertRdsToH5adFileOutputModel


class JobConvertRdsToH5adFileOutputModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobConvertRdsToH5adFileOutputModel
        fields = '__all__'
