import os
import subprocess

from celery import shared_task
from celery.exceptions import Reject, TaskError


@shared_task(bind=True)
def TreatmentVsControl(self, script_file, control_dataset_file, condition1_dataset_file, condition1_marker_genes_file, condition1_top25_markers_file, condition1_top10_genes_dotplot_file, condition2_marker_genes_file, condition2_top25_markers_file, condition2_top10_genes_dotplot_file, control_vs_conditions_markers_file, conditions_vs_control_markers_file, stdout_file, stderr_file, condition2_dataset_file=None):
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
            output_folder=os.path.dirname(conditions_vs_control_markers_file),
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
