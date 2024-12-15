from rest_framework import serializers

from ..models.JobInferenceModel import JobInferenceModel


class JobInferenceModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobInferenceModel
        fields = '__all__'
