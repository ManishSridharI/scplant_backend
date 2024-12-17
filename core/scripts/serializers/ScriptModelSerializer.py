from rest_framework import serializers

from ..models.ScriptModel import ScriptModel


class ScriptModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = ScriptModel
        fields = '__all__'
