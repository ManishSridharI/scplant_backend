from rest_framework import serializers

from ..models.H5adDatasetModel import H5adDatasetModel


class H5adDatasetModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = H5adDatasetModel
        fields = '__all__'
