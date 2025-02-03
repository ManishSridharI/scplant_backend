import os
from django.core.files.storage import FileSystemStorage


class OverwriteStorage(FileSystemStorage):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_available_name(self, name, max_length=None):
        if self.exists(name):
            self.delete(name)
        return super().get_available_name(name, max_length)

    def _save(self, name, content):
        try:
            os.makedirs(os.path.dirname(self.path(name)), exist_ok=True)
            with open(self.path(name), 'wb') as f:
                for chunk in content.chunks():
                    f.write(chunk)
        except Exception as e:
            print(e)
        return name
