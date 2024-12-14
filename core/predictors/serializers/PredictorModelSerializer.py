from rest_framework import serializers

from ..models.PredictorModel import PredictorModel


class PredictorModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = PredictorModel
        fields = '__all__'
