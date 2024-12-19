from rest_framework import serializers

from ..models.JobTreatmentVsControlModel import JobTreatmentVsControlModel


class JobTreatmentVsControlModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobTreatmentVsControlModel
        fields = '__all__'
