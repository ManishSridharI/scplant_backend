from rest_framework import serializers

from ..models.JobAnnotateAndPlotArgumentModel import JobAnnotateAndPlotArgumentModel


class JobAnnotateAndPlotArgumentModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobAnnotateAndPlotArgumentModel
        fields = '__all__'
