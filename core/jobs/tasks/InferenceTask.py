import subprocess

from celery import shared_task
from celery.exceptions import Reject, TaskError


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
