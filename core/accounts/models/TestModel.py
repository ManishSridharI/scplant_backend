import json
from django.db import models


class TestModel(models.Model):
    test = models.CharField(max_length=200, null=False, blank=False)
    test_description = models.TextField()

    class Meta:
        managed = False

    def __str__(self):
        output_string = json.dumps({
            "test": self.test,
            "test_description": self.test_description
        })
        return output_string

    @property
    def get_test(self):
        return f"{self.test}"

    @property
    def get_test_description(self):
        return f"{self.test_description}"
