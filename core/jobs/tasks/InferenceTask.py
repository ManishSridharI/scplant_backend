import subprocess

from celery import shared_task


@shared_task(bind=True)
def Inference(self, gene_number, dataset_file, predictor_file, log_file, prediction_file):
    try:
        print(gene_number)
        print(dataset_file)
        print(predictor_file)
        print(log_file)
        print(prediction_file)
        command = """
            pwd > uploads/current_directory.txt
        """
        result = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True
        )
        result = True
        return result
    except Exception as e:
        raise self.retry(exc=e, countdown=5, max_retries=2)
