import logging
from os import chdir
from os.path import join, basename
from celery import shared_task
from tempfile import mkdtemp
from shutil import copyfile, rmtree
from zipfile import ZipFile, is_zipfile
from prody import writePDB

from .models import Job, Result
from .combine_structures import combine_structures
from .knn_average import KNN_Averager
from .minimization import Minimizer


PYTHON = "/home/anaconda/bin/python"


@shared_task()
def run_netavg_calculation(job_id):
    job = Job.objects.get(pk=job_id)
    temp_dir = mkdtemp()
    chdir(temp_dir)
    trajectory_path = join(temp_dir, basename(job.trajectory.path))
    copyfile(job.trajectory.path, trajectory_path)
    logging.info( "%s" % trajectory_path)

    input_parse_success = True
    network_avg_success = True

    try:
        if trajectory_path.endswith('.zip') and is_zipfile(trajectory_path):
            # unzip, validate, group files for NetAvg calculation
            extract_to = join(temp_dir, 'extracted')
            with ZipFile(trajectory_path, 'r') as z:
                pdb_files = [f for f in z.namelist() if f.endswith('.pdb')]
                z.extractall(path=extract_to, members=pdb_files)
            # combine individual .pdb's into one multi-frame pdb file
            combined_pdb = join(temp_dir, 'combined.pdb')
            combine_structures(extract_to, combined_pdb)
            trajectory_file = combined_pdb

        elif trajectory_path.endswith('.pdb'):
            trajectory_file = trajectory_path
        
        else:
            job.status = Job.STATUS.error
            job.error_message = "Expected .zip or .pdb, got %s" % \
                                basename(trajectory_path)
            raise RuntimeError(job.error_message)

    except Exception as e:
        error_message = str(e)
        logging.error(str(e))
        job.status = Job.STATUS.error
        job.error_message = str(e)
        create_error_result_obj(job, error_message)
        input_parse_success = False


    if input_parse_success:
        try:
            #  =====================================
            #  = Compute network average structure =
            #  =====================================
            knn_input_path = trajectory_file
            knn_output_path = join(temp_dir, "avg.pdb")
            averager = KNN_Averager(knn_input_path)
            average_structure = averager.calc_average(job.knn)
            writePDB(knn_output_path, average_structure)

        except Exception as e:
            error_message = str(e)
            logging.error(str(e))
            job.status = Job.STATUS.error
            job.error_message = str(e)
            create_error_result_obj(job, error_message)
            network_avg_success = False

    #  ======================
    #  = Minimize structure =
    #  ======================
    if input_parse_success and network_avg_success:
        gmx_log_path = join(temp_dir, "gmx_log.txt")
        try:
            domin_input_path1 = knn_output_path
            domin_input_path2 = knn_input_path
            domin_output_path = \
                join(temp_dir, "netavg_%s" % basename(trajectory_file))
            m = Minimizer(domin_input_path1, domin_input_path2)
            with open(gmx_log_path, 'w') as f:
                minimized_protein = \
                    m.run_minimization(posres_force_const=1000., output_stream=f)
            writePDB(domin_output_path, minimized_protein)

            create_result_obj(job, domin_output_path)
            job.status = Job.STATUS.done

        except Exception as e:
            error_message = open(gmx_log_path, 'r').readlines()
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
