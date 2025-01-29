from rest_framework import serializers

from ..models.TenxfbcmDatasetModel import TenxfbcmDatasetModel


class TenxfbcmDatasetModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = TenxfbcmDatasetModel
        fields = '__all__'
