#!/usr/bin/env python

from glob import glob
import os
import sys
import subprocess

def checkruns(runname_str, out_sacct=True):

    # out_sacct can't be turned off if using the script currently; will change.

    rundirs = glob(runname_str)

    failed_dict = {}

    for rundir in rundirs:
        failed_massdirs = []
        jobids = []
        trackdirs = glob(os.path.join(rundir, '*M_dir'))
        for trackdir in trackdirs:
            # Checks for 'Fatal Error' string in the SLURM *.e file of the run's tracks:
            try:
                errfile = glob(os.path.join(trackdir, '*.e'))[0]
            except IndexError:
                print("{:s} is missing a SLURM .e (i.e., error) file.".format(trackdir))
                continue

            # Open .e file:
            with open(errfile, 'r') as ef:
                errlines = ef.readlines()
                for line in errlines:
                    if 'Fatal Error' in line:
                        print("{:s} encountered a fatal error.".format(trackdir))
                        failed_massdirs.append(trackdir)

            if out_sacct:
                # .o file is used to gather the SLURM job id; used to check job status.
                try:
                    outfile = glob(os.path.join(trackdir, '*.o'))[0]
                except IndexError:
                    print("{:s} is missing a SLURM .o (i.e., output) file.".format(trackdir))
                    continue

                with open(outfile, 'r') as of:
                    outlines = of.readlines()
                    for line in outlines:
                        if 'SLURM JOB ID:' in line:
                            jobids.append(line.split()[-1])
                            break

        if out_sacct:
            # Display sacct info for run's jobs:
            jobids = ','.join(jobids)
            print(subprocess.check_output(['sacct', '--jobs={:s}'.format(jobids), '--format=JobID,JobName,Partition,AllocCPUS,State,ExitCode,MaxRSS,Elapsed']))

        # Print number of failures and which masses will be re-sumitted:
        print("{:s} had {:d} tracks with fatal errors.".format(targetdir_of(rundir), len(failed_massdirs)))
        for failed_massdir in failed_massdirs:
            print("{:.2f} Msol will be resubmitted.".format(massof(failed_massdir.split('/')[-1])))

        # dictionary has the run directories as keys and its failed track subdirs as a list for the items.
        failed_dict[rundir] = failed_massdirs
        print("===========================================================")

    return failed_dict

def massof(trackname):

    # converts MIST track dir names into their resp. masses.

    mass_str = trackname.split('M')[0]

    return float(mass_str)/100.0

def targetdir_of(path):
  
    # returns the final directory name of a given path.

    return path.split('/')[-1]

if __name__ == '__main__':

    """
        Will check status of jobs for a MIST run & may also output the SLURM sacct info for that run's jobs. Also checks for certain
    fail states of runs and will re-submit the failed job upon discovery of a failure (although only checks for one fail state atm).

        Usage:
            
            >> ./check_jobs.py feh_p0.15*vvcrit*

        This would check runs with names following the above wildcard.

        The current fail state checked for is related to some SLURM scheduling/allocation failure for the job possibly? It involves
    a certain file being deemed inaccessible, leading the job to fail during the initial steps of the MESA run for certain tracks.


    NOTE: doesn't actually re-submit jobs; this functionality is currently commented out below.

    """


    # string for the run names to check; may contain wildcards. E.g. 'feh*TP'.
    runname_str = sys.argv[1]

    #prepend MIST grid path to runs:
    runname_str = os.path.join(os.environ['MIST_GRID_DIR'], 'MIST_v1.0', runname_str)

    # now check the runs for failed tracks:
    failed_dict = checkruns(runname_str)

    # Go through various run directories found...
    for rundir in failed_dict:
        # If there were failed tracks, re-submit them.
        failed_massdirs = failed_dict[rundir]
        print("___________________________________________________________")
        print("In {:s} re-submitting {:d} jobs.".format(targetdir_of(rundir), len(failed_massdirs)))

        # Re-submit the failed tracks:
        for failed_massdir in failed_massdirs:
            os.chdir(failed_massdir)
            print("re-submitting {:s}...".format(targetdir_of(failed_massdir)))
            slurmfile = glob(os.path.join(failed_massdir, '*M_run.sh'))[0]
            print(slurmfile)
            #os.system("sbatch " + slurmfile)
            #os.chdir(os.environ['MIST_CODE_DIR'])

