from rest_framework import serializers

from ..models.DatasetModel import DatasetModel


class DatasetModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = DatasetModel
        fields = '__all__'
