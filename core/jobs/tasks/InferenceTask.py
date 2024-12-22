import os
import datetime
import subprocess

from django.utils import timezone
from django.db import transaction

from celery import shared_task
from celery.exceptions import Reject, TaskError

from ..models.JobInferenceModel import JobInferenceModel
from ..models.JobInferenceArgumentModel import JobInferenceArgumentModel
from ..models.JobInferenceFileOutputModel import JobInferenceFileOutputModel


@shared_task(bind=True)
def RemoveInferenceTaskRecords(self):
    deadline = timezone.now() - datetime.timedelta(days=14, hours=0, minutes=0, seconds=0)
    job_inference_instance = JobInferenceModel.objects.filter(
        job_creation_timestamp__lt=deadline
    )
    job_celery_task_id_array = []
    if job_inference_instance.exists():
        for job_inference in job_inference_instance:
            if job_inference.job_celery_task_id not in job_celery_task_id_array:
                job_celery_task_id_array.append(job_inference.job_celery_task_id)
        with transaction.atomic():
            job_inference_instance.delete()
            if job_celery_task_id_array:
                if len(job_celery_task_id_array) > 0:
                    for job_celery_task_id in job_celery_task_id_array:
                        try:
                            job_inference_argument_instance = JobInferenceArgumentModel.objects.get(
                                job_celery_task_id=job_celery_task_id
                            )
                            job_inference_argument_instance.delete()
                        except Exception as e:
                            print(e)
                        try:
                            job_inference_file_output_instance = JobInferenceFileOutputModel.objects.get(
                                job_celery_task_id=job_celery_task_id
                            )
                            job_inference_file_output_instance.delete()
                        except Exception as e:
                            print(e)


@shared_task(bind=True)
def Inference(self, script_file, dataset_file, predictor_file, gene_number, log_file, prediction_file, stdout_file, stderr_file):
    try:
        if str(script_file).endswith(".py"):
            program = "python3"
        elif str(script_file).endswith(".R"):
            program = "Rscript"
        else:
            raise Reject(reason="Script type not supported", requeue=False)

        command = """
            {program} {script_file} \
            --data_path {dataset_file} \
            --model_path {predictor_file} \
            --gene_num {gene_number} \
            --log_file {log_file} \
            --prediction_file {prediction_file} > \
            {stdout_file} 2> \
            {stderr_file}
        """.format(
            program=program,
            script_file=script_file,
            dataset_file=dataset_file,
            predictor_file=predictor_file,
            gene_number=gene_number,
            log_file=log_file,
            prediction_file=prediction_file,
            stdout_file=stdout_file,
            stderr_file=stderr_file
        )

        command = command.replace("\n", " ")
        command = command.replace("\t", " ")

        completed_process_instance = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True
        )

        completed_process_instance.check_returncode()
        if completed_process_instance.returncode == 0:
            return True
        else:
            raise TaskError()
    except subprocess.CalledProcessError as e:
        raise self.retry(exc=e, countdown=5, max_retries=2)
    except Exception as e:
        raise self.retry(exc=e, countdown=5, max_retries=2)
