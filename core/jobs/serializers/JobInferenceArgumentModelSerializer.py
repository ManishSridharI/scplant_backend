from rest_framework import serializers

from ..models.JobInferenceArgumentModel import JobInferenceArgumentModel


class JobInferenceArgumentModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobInferenceArgumentModel
        fields = '__all__'
