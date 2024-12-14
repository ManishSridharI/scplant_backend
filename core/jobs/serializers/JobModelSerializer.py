from rest_framework import serializers

from ..models.JobModel import JobModel


class JobModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobModel
        fields = '__all__'
