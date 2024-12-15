from rest_framework import serializers

from ..models.JobAnnotateAndPlotFileOutputModel import JobAnnotateAndPlotFileOutputModel


class JobAnnotateAndPlotFileOutputModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobAnnotateAndPlotFileOutputModel
        fields = '__all__'
