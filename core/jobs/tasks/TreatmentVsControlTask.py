import os
import datetime
import subprocess

from django.utils import timezone
from django.db import transaction

from celery import shared_task
from celery.exceptions import Reject, TaskError

from ..models.JobTreatmentVsControlModel import JobTreatmentVsControlModel
from ..models.JobTreatmentVsControlFileOutputModel import JobTreatmentVsControlFileOutputModel


@shared_task(bind=True)
def RemoveTreatmentVsControlTaskRecords(self):
    deadline = timezone.now() - datetime.timedelta(days=14, hours=0, minutes=0, seconds=0)
    job_treatment_vs_control_instance = JobTreatmentVsControlModel.objects.filter(
        job_creation_timestamp__lt=deadline
    )
    job_celery_task_id_array = []
    if job_treatment_vs_control_instance.exists():
        for job_treatment_vs_control in job_treatment_vs_control_instance:
            if job_treatment_vs_control.job_celery_task_id not in job_celery_task_id_array:
                job_celery_task_id_array.append(job_treatment_vs_control.job_celery_task_id)
        with transaction.atomic():
            job_treatment_vs_control_instance.delete()
            if job_celery_task_id_array:
                if len(job_celery_task_id_array) > 0:
                    for job_celery_task_id in job_celery_task_id_array:
                        try:
                            job_treatment_vs_control_file_output_instance = JobTreatmentVsControlFileOutputModel.objects.get(
                                job_celery_task_id=job_celery_task_id
                            )
                            job_treatment_vs_control_file_output_instance.delete()
                        except Exception as e:
                            print(e)


@shared_task(bind=True)
def TreatmentVsControl(self, script_file, control_dataset_file, condition1_dataset_file, condition1_marker_genes_file, condition1_top25_markers_file, condition1_top10_genes_dotplot_file, stdout_file, stderr_file, condition2_dataset_file=None, condition2_marker_genes_file=None, condition2_top25_markers_file=None, condition2_top10_genes_dotplot_file=None, control_vs_conditions_markers_file=None, conditions_vs_control_markers_file=None):
    try:
        if str(script_file).endswith(".py"):
            program = "python3"
        elif str(script_file).endswith(".R"):
            program = "Rscript"
        else:
            raise Reject(reason="Script type not supported", requeue=False)

        command = """
            {program} {script_file} \
            --control_data_path {control_dataset_file} \
            --condition1_data_path {condition1_dataset_file} \
        """.format(
            program=program,
            script_file=script_file,
            control_dataset_file=control_dataset_file,
            condition1_dataset_file=condition1_dataset_file
        )

        if condition2_dataset_file:
            command += """
                --condition2_data_path {condition2_dataset_file} \
            """.format(
                condition2_dataset_file=condition2_dataset_file
            )

        command += """
            --output_folder {output_folder} > \
            {stdout_file} 2> \
            {stderr_file}
        """.format(
            output_folder=os.path.dirname(condition1_marker_genes_file),
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
