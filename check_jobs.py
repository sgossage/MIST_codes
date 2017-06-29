#!/usr/bin/env python

from glob import glob
import os
import sys

def checkruns(runname_str):

    rundirs = glob(runname_str)

    failed_dict = {}

    for rundir in rundirs:
        failed_massdirs = []
        trackdirs = glob(os.path.join(rundir, '*M_dir'))
        for trackdir in trackdirs:
            # Checks for 'Fatal Error' string in the SLURM *.e file of the run's tracks:
            try:
                errfile = glob(os.path.join(trackdir, '*.e'))[0]
            except IndexError:
                print("{:s} is missing a SLURM .e (i.e., error) file.".format(trackdir))
                continue

            with open(errfile, 'r') as ef:
                errlines = ef.readlines()
                for line in errlines:
                    if 'Fatal Error' in line:
                        print("{:s} encountered a fatal error.".format(trackdir))
                        #failed_masses.append(massof(trackdir.split('/')[-1]))
                        failed_massdirs.append(trackdir)
        print("{:s} had {:d} tracks with fatal errors.".format(rundir, len(failed_massdirs)))
        for failed_massdir in failed_massdirs:
            print("{:.2f} Msol will be resubmitted.".format(massof(failed_massdir.split('/')[-1])))

        # dictionary has the run directories as keys and its failed track subdirs as a list for the items.
        failed_dict[rundir] = failed_massdirs

    return failed_dict

def massof(trackname):

    # converts MIST track dir names into their resp. masses.

    mass_str = trackname.split('M')[0]

    return float(mass_str)/100.0

def targetdir_of(path):

    return path.split('/')[-1]

if __name__ == '__main__':

    # string for the run names to check; may contain wildcards. E.g. 'feh*TP'.
    runname_str = sys.argv[1]

    #prepend MIST grid path to runs:
    runname_str = os.path.join(os.environ['MIST_GRID_DIR'], 'MIST_v1.0', runname_str)

    # now check the runs for failed tracks:
    failed_dict = checkruns(runname_str)

    for rundir in failed_dict:
        failed_massdirs = failed_dict[rundir]
        print("___________________________________________________________")
        print("In {:s} re-submitting {:d} jobs.".format(targetdir_of(rundir), len(failed_massdirs)))
        for failed_massdir in failed_massdirs:
            os.chdir(failed_massdir)
            print("re-submitting {:s}...".format(targetdir_of(failed_massdir)))
            slurmfile = glob(os.path.join(failed_massdir, '*M_run.sh'))[0]
            print(slurmfile)
            #os.system("sbatch " + slurmfile)
            #os.chdir(os.environ['MIST_CODE_DIR'])

