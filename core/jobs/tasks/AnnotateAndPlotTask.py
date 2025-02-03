import os
import datetime
import subprocess

from django.utils import timezone
from django.db import transaction

from celery import shared_task
from celery.exceptions import Reject, TaskError

from ..models.JobAnnotateAndPlotModel import JobAnnotateAndPlotModel
from ..models.JobAnnotateAndPlotFileOutputModel import JobAnnotateAndPlotFileOutputModel


@shared_task(bind=True)
def RemoveAnnotateAndPlotTaskRecords(self):
    deadline = timezone.now() - datetime.timedelta(days=14, hours=0, minutes=0, seconds=0)
    job_annotate_and_plot_instance = JobAnnotateAndPlotModel.objects.filter(
        job_creation_timestamp__lt=deadline
    )
    job_celery_task_id_array = []
    if job_annotate_and_plot_instance.exists():
        for job_annotate_and_plot in job_annotate_and_plot_instance:
            if job_annotate_and_plot.job_celery_task_id not in job_celery_task_id_array:
                job_celery_task_id_array.append(job_annotate_and_plot.job_celery_task_id)
        with transaction.atomic():
            job_annotate_and_plot_instance.delete()
            if job_celery_task_id_array:
                if len(job_celery_task_id_array) > 0:
                    for job_celery_task_id in job_celery_task_id_array:
                        try:
                            job_annotate_and_plot_file_output_instance = JobAnnotateAndPlotFileOutputModel.objects.get(
                                job_celery_task_id=job_celery_task_id
                            )
                            job_annotate_and_plot_file_output_instance.delete()
                        except Exception as e:
                            print(e)


@shared_task(bind=True)
def AnnotateAndPlot(self, script_file, h5ad_dataset_file, tenxfbcm_dataset_folder, pred_dataset_file, predictor_file, data_type, folder_path, output_prediction_file, output_log_file, output_stats_csv_file, stdout_file, stderr_file, output_pred_dataset_file):
    try:
        if str(script_file).endswith(".py"):
            program = "python3"
        elif str(script_file).endswith(".R"):
            program = "Rscript"
        else:
            raise Reject(reason="Script type not supported", requeue=False)

        command = " "
        command = command + "rm -rf {output_stats_csv_file}; ".format(
            output_stats_csv_file=output_stats_csv_file
        )
        if pred_dataset_file:
            command = command + "cp -r {input_file} {output_file}; ".format(
                input_file=pred_dataset_file,
                output_file=output_prediction_file
            )
        else:
            command = command + "rm -rf {output_prediction_file}; ".format(
                output_prediction_file=output_prediction_file
            )
        command = command + "{program} {script_file} ".format(program=program, script_file=script_file)
        command = command + "--model_path {predictor_file} ".format(predictor_file=predictor_file)
        if h5ad_dataset_file:
            command = command + "--data_path {data_path} ".format(data_path=h5ad_dataset_file)
        if tenxfbcm_dataset_folder:
            command = command + "--data_path {data_path} ".format(data_path=tenxfbcm_dataset_folder)
        command = command + "--data_type {data_type} ".format(data_type=data_type)
        command = command + "--output_folder {output_folder} ".format(output_folder=folder_path)
        command = command + "--log {output_log_file} ".format(output_log_file=output_log_file)
        command = command + "--prediction_file {prediction_file} ".format(prediction_file=os.path.basename(output_prediction_file))
        command = command + "--stats_file {stats_csv_file} ".format(stats_csv_file=os.path.basename(output_stats_csv_file))
        command = command + "> {stdout_file} 2> {stderr_file}; ".format(stdout_file=stdout_file, stderr_file=stderr_file)
        command = command + "cp -r {prediction_file} {output_pred_dataset_file}; ".format(prediction_file=output_prediction_file, output_pred_dataset_file=output_pred_dataset_file)

        command = command.replace("\n", " ")
        command = command.replace("\t", " ")

        print(command)

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
