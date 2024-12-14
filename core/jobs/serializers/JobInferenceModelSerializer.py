from rest_framework import serializers

from django.core.files.base import ContentFile, File

from ..models.JobInferenceModel import JobInferenceModel


class JobInferenceModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobInferenceModel
        fields = '__all__'
