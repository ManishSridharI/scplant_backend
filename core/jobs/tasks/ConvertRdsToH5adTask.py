import os
import datetime
import subprocess

from django.utils import timezone
from django.db import transaction

from celery import shared_task
from celery.exceptions import Reject, TaskError

from ..models.JobConvertRdsToH5adModel import JobConvertRdsToH5adModel
from ..models.JobConvertRdsToH5adFileOutputModel import JobConvertRdsToH5adFileOutputModel


@shared_task(bind=True)
def RemoveConvertRdsToH5adTaskRecords(self):
    deadline = timezone.now() - datetime.timedelta(days=30, hours=0, minutes=0, seconds=0)
    job_convert_rds_to_h5ad_instance = JobConvertRdsToH5adModel.objects.filter(
        job_creation_timestamp__lt=deadline
    )
    job_celery_task_id_array = []
    if job_convert_rds_to_h5ad_instance.exists():
        for job_convert_rds_to_h5ad in job_convert_rds_to_h5ad_instance:
            if job_convert_rds_to_h5ad.job_celery_task_id not in job_celery_task_id_array:
                job_celery_task_id_array.append(job_convert_rds_to_h5ad.job_celery_task_id)
        with transaction.atomic():
            job_convert_rds_to_h5ad_instance.delete()
            if job_celery_task_id_array:
                if len(job_celery_task_id_array) > 0:
                    for job_celery_task_id in job_celery_task_id_array:
                        try:
                            job_convert_rds_to_h5ad_file_output_instance = JobConvertRdsToH5adFileOutputModel.objects.get(
                                job_celery_task_id=job_celery_task_id
                            )
                            job_convert_rds_to_h5ad_file_output_instance.delete()
                        except Exception as e:
                            print(e)


@shared_task(bind=True)
def ConvertRdsToH5ad(self, script_file, rds_dataset, output_file, stdout_file, stderr_file, h5ad_dataset):
    try:
        if str(script_file).endswith(".py"):
            program = "python3"
        elif str(script_file).endswith(".R"):
            program = "Rscript"
        else:
            raise Reject(reason="Script type not supported", requeue=False)

        command = """
            {program} {script_file} \
            --input {rds_dataset} \
            --output {output_file}  > \
            {stdout_file} 2> \
            {stderr_file};
            cp -rf {output_file} {h5ad_dataset};
        """.format(
            program=program,
            script_file=script_file,
            rds_dataset=rds_dataset,
            output_file=output_file,
            stdout_file=stdout_file,
            stderr_file=stderr_file,
            h5ad_dataset=h5ad_dataset
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
