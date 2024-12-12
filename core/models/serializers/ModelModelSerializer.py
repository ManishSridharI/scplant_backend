from rest_framework import serializers

from ..models.ModelModel import ModelModel


class ModelModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = ModelModel
        fields = '__all__'
