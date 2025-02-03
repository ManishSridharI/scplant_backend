from rest_framework import serializers

from ..models.JobTreatmentVsControlFileOutputModel import JobTreatmentVsControlFileOutputModel


class JobTreatmentVsControlFileOutputModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobTreatmentVsControlFileOutputModel
        fields = '__all__'
