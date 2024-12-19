from rest_framework import serializers

from ..models.JobCompareCellTypeDistFileOutputModel import JobCompareCellTypeDistFileOutputModel


class JobCompareCellTypeDistFileOutputModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobCompareCellTypeDistFileOutputModel
        fields = '__all__'
