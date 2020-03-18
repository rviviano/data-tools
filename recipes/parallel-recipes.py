from __future__ import print_function
import os, sys, time, traceback
import multiprocessing as mp
from multiprocessing import Process, Lock


def target_function():
    """ Function to process data. Called by run_parallel """
    pass


def run_embarass_parallel(data_to_process, cores=24):
    """ General code for embarassingly parallel python multiprocessing. E.g.,
        code to apply the same processing to multiple datasets with the same
        structure. Spawns a new process for each individual data to process and
        processes each unique dataset separately."""

    # Process lock for stdout access
    print_lock = Lock()

    # Dictionary to keep track of running processes
    ps = {} 

    for x in data_to_process:
        # May need to apply preprocessing prior to argument determination here
        try:
            # Specify arguments to target function
            arguments = (x, print_lock)
            # Name the process
            name = ()
            # Specify and start the process
            p = Process(name=name, target=target_function, args=arguments)
            p.start()
            ps[name] = p
        except KeyboardInterrupt:
            raise
        except:
            traceback.print_exc(file=sys.stdout)
            pass

        # Check if any process has completed, joins a process if it has and then 
        # exits this while loop to allow a new process to start
        while len(ps.keys()) >= (cores):  
            for j in ps.keys():
                if not ps[j].is_alive():
                    ps[j].join()
                    del ps[j]
            time.sleep(5)