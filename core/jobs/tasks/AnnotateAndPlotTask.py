import os
import datetime
import subprocess

from django.utils import timezone
from django.db import transaction

from celery import shared_task
from celery.exceptions import Reject, TaskError

from ..models.JobAnnotateAndPlotModel import JobAnnotateAndPlotModel
from ..models.JobAnnotateAndPlotArgumentModel import JobAnnotateAndPlotArgumentModel
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
                            job_annotate_and_plot_argument_instance = JobAnnotateAndPlotArgumentModel.objects.get(
                                job_celery_task_id=job_celery_task_id
                            )
                            job_annotate_and_plot_argument_instance.delete()
                        except Exception as e:
                            print(e)
                        try:
                            job_annotate_and_plot_file_output_instance = JobAnnotateAndPlotFileOutputModel.objects.get(
                                job_celery_task_id=job_celery_task_id
                            )
                            job_annotate_and_plot_file_output_instance.delete()
                        except Exception as e:
                            print(e)


@shared_task(bind=True)
def AnnotateAndPlot(self, script_file, dataset_file, predictor_file, gene_number, top25_markers_file, marker_genes_file, output_with_celltype_file, log_file, prediction_file, annotate_tsne_file, annotate_umap_file, top3_genes_dotplot_file, stdout_file, stderr_file):
    try:
        if str(script_file).endswith(".py"):
            program = "python3"
        elif str(script_file).endswith(".R"):
            program = "Rscript"
        else:
            raise Reject(reason="Script type not supported", requeue=False)

        command = """
            rm -rf {prediction_file};
            {program} {script_file} \
            --data_path {dataset_file} \
            --model_path {predictor_file} \
            --gene_num {gene_number} \
            --log_file {log_file} \
            --output_folder {output_folder} > \
            {stdout_file} 2> \
            {stderr_file}
        """.format(
            prediction_file=prediction_file,
            program=program,
            script_file=script_file,
            dataset_file=dataset_file,
            predictor_file=predictor_file,
            gene_number=gene_number,
            log_file=log_file,
            output_folder=os.path.dirname(prediction_file),
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
