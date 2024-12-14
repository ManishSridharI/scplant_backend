import subprocess

from celery import shared_task


@shared_task(bind=True)
def Add(self, x, y):
    try:
        result = x + y
        print(result)
        return result
    except Exception as e:
        raise self.retry(exc=e, countdown=5, max_retries=2)


@shared_task(bind=True)
def WriteDate(self):
    try:
        command = """
            date > uploads/current_date.txt
        """
        result = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True
        )
    except Exception as e:
        raise self.retry(exc=e, countdown=5, max_retries=2)
