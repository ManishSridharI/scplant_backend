import subprocess

from celery import shared_task
from celery.signals import before_task_publish, after_task_publish, task_prerun, task_received, task_success, task_failure, task_internal_error, task_revoked, task_unknown, task_rejected
from celery.result import AsyncResult

from ..models.JobCompareCellTypeDistModel import JobCompareCellTypeDistModel

from ..serializers.JobCompareCellTypeDistModelSerializer import JobCompareCellTypeDistModelSerializer

from ..tasks.CompareCellTypeDistTask import CompareCellTypeDist


def update_job_compare_cell_type_dist_celery_task_status_and_result(job_celery_task_id, job_celery_task_status, job_celery_task_result=None):
    try:
        job_compare_cell_type_dist_model_instance = JobCompareCellTypeDistModel.objects.get(
            job_celery_task_id=job_celery_task_id
        )
        if job_compare_cell_type_dist_model_instance:
            if job_celery_task_status and job_celery_task_result:
                job_compare_cell_type_dist_model_serializer = JobCompareCellTypeDistModelSerializer(
                    instance=job_compare_cell_type_dist_model_instance,
                    data={
                        'job_celery_task_status': job_celery_task_status,
                        'job_celery_task_result': job_celery_task_result
                    },
                    partial=True
                )
            else:
                job_compare_cell_type_dist_model_serializer = JobCompareCellTypeDistModelSerializer(
                    instance=job_compare_cell_type_dist_model_instance,
                    data={
                        'job_celery_task_status': job_celery_task_status
                    },
                    partial=True
                )
            if job_compare_cell_type_dist_model_serializer.is_valid():
                try:
                    job_compare_cell_type_dist_model_serializer_instance = job_compare_cell_type_dist_model_serializer.save()
                except Exception as e:
                    print(e)
            else:
                print(job_compare_cell_type_dist_model_serializer.errors)
    except Exception as e:
        print(e)


@task_received.connect(sender=CompareCellTypeDist)
def update_on_task_received(sender=None, result=None, **kwargs):
    update_job_compare_cell_type_dist_celery_task_status_and_result(
        job_celery_task_id=str(sender.request.id),
        job_celery_task_status="PENDING"
    )


@task_prerun.connect(sender=CompareCellTypeDist)
def update_on_task_prerun(sender=None, result=None, **kwargs):
    update_job_compare_cell_type_dist_celery_task_status_and_result(
        job_celery_task_id=str(sender.request.id),
        job_celery_task_status="RUNNING"
    )


@task_success.connect(sender=CompareCellTypeDist)
def update_on_task_success(sender=None, result=None, **kwargs):
    async_result_instance = AsyncResult(str(sender.request.id))
    update_job_compare_cell_type_dist_celery_task_status_and_result(
        job_celery_task_id=str(sender.request.id),
        job_celery_task_status=str(async_result_instance.status),
        job_celery_task_result=str(result)
    )


@task_failure.connect(sender=CompareCellTypeDist)
def update_on_task_failure(sender=None, exception=None, **kwargs):
    async_result_instance = AsyncResult(str(sender.request.id))
    update_job_compare_cell_type_dist_celery_task_status_and_result(
        job_celery_task_id=str(sender.request.id),
        job_celery_task_status=str(async_result_instance.status),
        job_celery_task_result=str(exception)
    )


@task_internal_error.connect(sender=CompareCellTypeDist)
def update_on_task_internal_error(sender=None, result=None, **kwargs):
    update_job_compare_cell_type_dist_celery_task_status_and_result(
        job_celery_task_id=str(sender.request.id),
        job_celery_task_status="INTERNAL ERROR"
    )


@task_revoked.connect(sender=CompareCellTypeDist)
def update_on_task_revoked(sender=None, result=None, **kwargs):
    update_job_compare_cell_type_dist_celery_task_status_and_result(
        job_celery_task_id=str(sender.request.id),
        job_celery_task_status="REVOKED"
    )


@task_unknown.connect(sender=CompareCellTypeDist)
def update_on_task_unknown(sender=None, result=None, **kwargs):
    update_job_compare_cell_type_dist_celery_task_status_and_result(
        job_celery_task_id=str(sender.request.id),
        job_celery_task_status="UNKNOWN"
    )


@task_rejected.connect(sender=CompareCellTypeDist)
def update_on_task_rejected(sender=None, result=None, **kwargs):
    update_job_compare_cell_type_dist_celery_task_status_and_result(
        job_celery_task_id=str(sender.request.id),
        job_celery_task_status="REJECTED"
    )
