import logging
import subprocess
from os import chdir
from os.path import join, basename
from celery import shared_task
from tempfile import mkdtemp
from shutil import copyfile, rmtree
from zipfile import ZipFile, is_zipfile
from glob import glob

from .models import Job, Result

@shared_task()
def run_netavg_calculation(job_id):
    job = Job.objects.get(pk=job_id)
    temp_dir = mkdtemp()
    chdir(temp_dir)
    trajectory_path = join(temp_dir, basename(job.trajectory.path))
    copyfile(job.trajectory.path, trajectory_path)
    logging.info( "%s" % trajectory_path)

    if trajectory_path.endswith('.zip') and is_zipfile(trajectory_path):
        # unzip, validate, group files for NetAvg calculation
        extract_to = join(temp_dir, 'extracted')
        with ZipFile(trajectory_path, 'r') as z:
            pdb_files = [f for f in z.namelist() if f.endswith('.pdb')]
            z.extractall(path=extract_to, members=pdb_files)
        pdb_files = glob( join(extract_to, '*.pdb') )
        # combine individual .pdb's into one multi-frame pdb file
        trajectory_file = None

    elif trajectory_path.endswith('.pdb'):
        trajectory_file = trajectory_path
    
    else:
        job.status = Job.STATUS.error
        job.error_message = "Expected .zip or .pdb, got %s" % \
                            basename(trajectory_path)

    knn_output_str = ''
    domin_output_str = ''
    try:
        python_cmd = "/home/anaconda/bin/python"
        knn_cmd = "/home/NetAvg/knn_average.py"
        knn_option = str(job.knn)
        knn_input_path = trajectory_file
        knn_output_path = join(temp_dir, "avg.pdb")
        arg_list = [python_cmd, knn_cmd, '--knn', knn_option, knn_input_path,
                    knn_output_path]
        knn_output_str = \
            subprocess.check_output(arg_list, stderr=subprocess.STDOUT)

        domin_cmd = "/home/NetAvg/do_minimization.py"
        domin_input_path1 = knn_output_path
        domin_input_path2 = knn_input_path
        domin_output_path = join(temp_dir, "netavg_%s" % basename(trajectory_file))
        arg_list2 = [python_cmd, domin_cmd, domin_input_path1, domin_input_path2,
                     domin_output_path]
        domin_output_str = \
            subprocess.check_output(arg_list2, stderr=subprocess.STDOUT)

        create_result_obj(job, domin_output_path)
        # after all comparison files have been processed
        job.status = Job.STATUS.done

    except Exception as e:
        error_message = knn_output_str + '\n' + domin_output_str
        logging.error(str(e))
        job.status = Job.STATUS.error
        job.error_message = str(e)
        create_error_result_obj(job, error_message)

    job.save()
    rmtree(temp_dir)
    job.trajectory.delete()
    return

def create_result_obj(job, output_path):
    print job, output_path
    with open(output_path, 'r') as f:
        output_str = f.read()
    r = Result.objects.create(job=job, output_pdb=output_str)
    return

def create_error_result_obj(job, error_message):
    r = Result.objects.create(job=job, error_message=error_message)
    return
