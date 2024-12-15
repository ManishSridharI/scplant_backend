import subprocess

from celery import shared_task


@shared_task(bind=True)
def Inference(self, gene_number, dataset_file, predictor_file, log_file, prediction_file, stdout_file, stderr_file):
    try:
        print(gene_number)
        print(dataset_file)
        print(predictor_file)
        print(log_file)
        print(prediction_file)
        print(stdout_file)
        print(stderr_file)
        command = """
            python3 ../model_codebase/inference.py \
            --gene_num {gene_number} \
            --data_path {dataset_file} \
            --model_path {predictor_file} \
            --log_file {log_file} \
            --prediction_file {prediction_file} > \
            {stdout_file} 2> \
            {stderr_file}
        """.format(
            gene_number=gene_number,
            dataset_file=dataset_file,
            predictor_file=predictor_file,
            log_file=log_file,
            prediction_file=prediction_file,
            stdout_file=stdout_file,
            stderr_file=stderr_file
        )
        print(command)
        completed_process_instance = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True
        )
        completed_process_instance.check_returncode()
        return True
    except subprocess.CalledProcessError as e:
        raise self.retry(exc=e, countdown=5, max_retries=2)
    except Exception as e:
        raise self.retry(exc=e, countdown=5, max_retries=2)
