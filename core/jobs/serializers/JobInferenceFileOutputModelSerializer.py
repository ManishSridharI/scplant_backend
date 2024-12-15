from rest_framework import serializers

from ..models.JobInferenceFileOutputModel import JobInferenceFileOutputModel


class JobInferenceFileOutputModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobInferenceFileOutputModel
        fields = '__all__'
