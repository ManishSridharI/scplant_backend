from rest_framework import serializers

from ..models.JobCompareCellTypeDistModel import JobCompareCellTypeDistModel


class JobCompareCellTypeDistModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobCompareCellTypeDistModel
        fields = '__all__'
