from rest_framework import serializers

from ..models.OrganismModel import OrganismModel


class OrganismModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = OrganismModel
        fields = '__all__'
