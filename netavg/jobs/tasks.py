import logging
import subprocess
from os import chdir
from os.path import join, basename
from celery import shared_task
from tempfile import mkdtemp
from shutil import copyfile, rmtree
from zipfile import ZipFile, is_zipfile

from .models import Job, Result


PYTHON = "/home/anaconda/bin/python"


@shared_task()
def run_netavg_calculation(job_id):
    job = Job.objects.get(pk=job_id)
    temp_dir = mkdtemp()
    chdir(temp_dir)
    trajectory_path = join(temp_dir, basename(job.trajectory.path))
    copyfile(job.trajectory.path, trajectory_path)
    logging.info( "%s" % trajectory_path)

    combine_output_str = ''
    knn_output_str = ''
    domin_output_str = ''

    try:

        if trajectory_path.endswith('.zip') and is_zipfile(trajectory_path):
            # unzip, validate, group files for NetAvg calculation
            extract_to = join(temp_dir, 'extracted')
            with ZipFile(trajectory_path, 'r') as z:
                pdb_files = [f for f in z.namelist() if f.endswith('.pdb')]
                z.extractall(path=extract_to, members=pdb_files)
            # combine individual .pdb's into one multi-frame pdb file
            combine_pdb_cmd = "/home/NetAvg/combine_structures.py"
            combined_pdb = join(temp_dir, 'combined.pdb')
            arg_list = [PYTHON, combine_pdb_cmd, extract_to, combined_pdb]
            combine_output_str = \
                subprocess.check_output(arg_list, stderr=subprocess.STDOUT)
            trajectory_file = combined_pdb

        elif trajectory_path.endswith('.pdb'):
            trajectory_file = trajectory_path
        
        else:
            job.status = Job.STATUS.error
            job.error_message = "Expected .zip or .pdb, got %s" % \
                                basename(trajectory_path)

        knn_cmd = "/home/NetAvg/knn_average.py"
        knn_option = str(job.knn)
        knn_input_path = trajectory_file
        knn_output_path = join(temp_dir, "avg.pdb")
        arg_list = [PYTHON, knn_cmd, '--knn', knn_option, knn_input_path,
                    knn_output_path]
        knn_output_str = \
            subprocess.check_output(arg_list, stderr=subprocess.STDOUT)

        domin_cmd = "/home/NetAvg/do_minimization.py"
        domin_input_path1 = knn_output_path
        domin_input_path2 = knn_input_path
        domin_output_path = join(temp_dir, "netavg_%s" % basename(trajectory_file))
        arg_list2 = [PYTHON, domin_cmd, domin_input_path1, domin_input_path2,
                     domin_output_path]
        domin_output_str = \
            subprocess.check_output(arg_list2, stderr=subprocess.STDOUT)

        create_result_obj(job, domin_output_path)
        # after all comparison files have been processed
        job.status = Job.STATUS.done

    except Exception as e:
        error_message = combine_output_str + '\n' + knn_output_str + '\n' + domin_output_str
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
