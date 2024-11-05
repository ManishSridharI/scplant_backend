from rest_framework import serializers

from ..models.CustomUserModel import CustomUserModel


class CustomUserModelSerializer(serializers.ModelSerializer):
    class Meta:
        model = CustomUserModel
        fields = '__all__'
        extra_kwargs = {
            'password': {'write_only': True}
        }
