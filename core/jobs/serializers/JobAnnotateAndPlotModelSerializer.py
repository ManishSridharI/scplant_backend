from rest_framework import serializers

from ..models.JobAnnotateAndPlotModel import JobAnnotateAndPlotModel


class JobAnnotateAndPlotModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobAnnotateAndPlotModel
        fields = '__all__'
