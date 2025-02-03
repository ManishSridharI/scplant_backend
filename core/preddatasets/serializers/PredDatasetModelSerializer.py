from rest_framework import serializers

from ..models.PredDatasetModel import PredDatasetModel


class PredDatasetModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = PredDatasetModel
        fields = '__all__'
