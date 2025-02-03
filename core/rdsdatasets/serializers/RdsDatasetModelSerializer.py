from rest_framework import serializers

from ..models.RdsDatasetModel import RdsDatasetModel


class RdsDatasetModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = RdsDatasetModel
        fields = '__all__'
