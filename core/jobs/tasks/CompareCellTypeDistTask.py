import os
import datetime
import subprocess

from django.utils import timezone
from django.db import transaction

from celery import shared_task
from celery.exceptions import Reject, TaskError

from ..models.JobCompareCellTypeDistModel import JobCompareCellTypeDistModel
from ..models.JobCompareCellTypeDistFileOutputModel import JobCompareCellTypeDistFileOutputModel


@shared_task(bind=True)
def RemoveCompareCellTypeDistTaskRecords(self):
    deadline = timezone.now() - datetime.timedelta(days=14, hours=0, minutes=0, seconds=0)
    job_compare_cell_type_dist_instance = JobCompareCellTypeDistModel.objects.filter(
        job_creation_timestamp__lt=deadline
    )
    job_celery_task_id_array = []
    if job_compare_cell_type_dist_instance.exists():
        for job_compare_cell_type_dist in job_compare_cell_type_dist_instance:
            if job_compare_cell_type_dist.job_celery_task_id not in job_celery_task_id_array:
                job_celery_task_id_array.append(job_compare_cell_type_dist.job_celery_task_id)
        with transaction.atomic():
            job_compare_cell_type_dist_instance.delete()
            if job_celery_task_id_array:
                if len(job_celery_task_id_array) > 0:
                    for job_celery_task_id in job_celery_task_id_array:
                        try:
                            job_compare_cell_type_dist_file_output_instance = JobCompareCellTypeDistFileOutputModel.objects.get(
                                job_celery_task_id=job_celery_task_id
                            )
                            job_compare_cell_type_dist_file_output_instance.delete()
                        except Exception as e:
                            print(e)


@shared_task(bind=True)
def CompareCellTypeDist(self, script_file, control_prediction_file, condition1_prediction_file, output_file, stdout_file, stderr_file, condition2_prediction_file=None):
    try:
        if str(script_file).endswith(".py"):
            program = "python3"
        elif str(script_file).endswith(".R"):
            program = "Rscript"
        else:
            raise Reject(reason="Script type not supported", requeue=False)

        command = """
            {program} {script_file} \
            --control_pred_file {control_prediction_file} \
            --condition1_pred_file {condition1_prediction_file} \
        """.format(
            program=program,
            script_file=script_file,
            control_prediction_file=control_prediction_file,
            condition1_prediction_file=condition1_prediction_file
        )

        if condition2_prediction_file:
            command += """
                --condition2_pred_file {condition2_prediction_file} \
            """.format(
                condition2_prediction_file=condition2_prediction_file
            )

        command += """
            --output_folder {output_folder} > \
            {stdout_file} 2> \
            {stderr_file}
        """.format(
            output_folder=os.path.dirname(output_file),
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
