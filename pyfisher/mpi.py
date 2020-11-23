from __future__ import print_function
import numpy as np
import os,sys,time

"""
Copied from orphics.mpi
"""

try:
    disable_mpi_env = os.environ['DISABLE_MPI']
    disable_mpi = True if disable_mpi_env.lower().strip() == "true" else False
except:
    disable_mpi = False

"""
Use the below cleanup stuff only for intel-mpi!
If you use it on openmpi, you will have no traceback for errors
causing hours of endless confusion and frustration! - Sincerely, past frustrated Mat
"""
# From Sigurd's enlib.mpi:
# Uncaught exceptions don't cause mpi to abort. This can lead to thousands of
# wasted CPU hours
# def cleanup(type, value, traceback):
#     sys.__excepthook__(type, value, traceback)
#     MPI.COMM_WORLD.Abort(1)
# sys.excepthook = cleanup


class fakeMpiComm:
    """
    A Simple Fake MPI implementation
    """
    def __init__(self):
        pass
    def Get_rank(self):
        return 0
    def Get_size(self):
        return 1
    def Barrier(self):
        pass
    def Abort(self,dummy):
        pass




try:
    if disable_mpi: raise
    from mpi4py import MPI
except:

    if not(disable_mpi): print("WARNING: mpi4py could not be loaded. Falling back to fake MPI. This means that if you submitted multiple processes, they will all be assigned the same rank of 0, and they are potentially doing the same thing.")
    
    class template:
        pass

    MPI = template()
    MPI.COMM_WORLD = fakeMpiComm()


def mpi_distribute(num_tasks,avail_cores,allow_empty=False):
    # copied to mapsims.convert_noise_templates
    if not(allow_empty): assert avail_cores<=num_tasks
    min_each, rem = divmod(num_tasks,avail_cores)
    num_each = np.array([min_each]*avail_cores) # first distribute equally
    if rem>0: num_each[-rem:] += 1  # add the remainder to the last set of cores (so that rank 0 never gets extra jobs)

    task_range = list(range(num_tasks)) # the full range of tasks
    cumul = np.cumsum(num_each).tolist() # the end indices for each task
    task_dist = [task_range[x:y] for x,y in zip([0]+cumul[:-1],cumul)] # a list containing the tasks for each core
    assert sum(num_each)==num_tasks
    assert len(num_each)==avail_cores
    assert len(task_dist)==avail_cores
    return num_each,task_dist
    


def distribute(njobs,verbose=True,**kwargs):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    numcores = comm.Get_size()
    num_each,each_tasks = mpi_distribute(njobs,numcores,**kwargs)
    if rank==0: print ("At most ", max(num_each) , " tasks...")
    my_tasks = each_tasks[rank]
    return comm,rank,my_tasks
