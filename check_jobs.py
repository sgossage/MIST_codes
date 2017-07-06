#!/usr/bin/env python

from glob import glob
import os
import sys
import subprocess
import shutil
from distutils.dir_util import copy_tree

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
        print("{:s} had {:d} tracks with fatal errors.".format(target_of(rundir), len(failed_massdirs)))
        for failed_massdir in failed_massdirs:
            print("{:.2f} Msol will be resubmitted.".format(massof(failed_massdir.split('/')[-1])))

        # dictionary has the run directories as keys and its failed track subdirs as a list for the items.
        failed_dict[rundir] = failed_massdirs
        print("===========================================================")

    return failed_dict

def getphotos(runname_str, TPextend=True):

    # out_sacct can't be turned off if using the script currently; will change.

    rundirs = glob(runname_str)

    res_dict = {}

    for rundir in rundirs:
        badtracks = []
        trackdirs = glob(os.path.join(rundir, '*M*dir'))
        for trackdir in trackdirs:
            # Checks for the final photo of a track:

            # .o file is used to recover the final photo of the MESA run.
            try:
                outfile = glob(os.path.join(trackdir, '*.o'))[0]
            except IndexError:
                print("{:s} is missing a SLURM .o (i.e., output) file.".format(trackdir))
                continue

            with open(outfile, 'r') as of:
                outlines = of.readlines()
                for i, line in enumerate(outlines):
                    if 'save photos/x' in line:
                        if TPextend == True & ('termination code: HB_limit' in outlines[i+1]) & ( not glob(os.path.join(trackdir, '*re.*'))):
                            # check if extension to TP is necc.
                            lastphoto = line.split('photos/')[-1].split(' ')[0]
                            #photos.append((target_of(trackdir), lastphoto))
                            badtracks.append(trackdir)
                            # Make a SLURM submission script to resume the job.
                            # first copying the XXXXXM_run.sh script to make a new script that will resume the run:
                            runscript_path = os.path.join(trackdir, target_of(trackdir).replace('dir', 'run.sh'))
                            resscript_path = os.path.join(trackdir, target_of(runscript_path).replace('run', 'resume'))
                            # modify the 'resume' script:
                            with open(resscript_path, 'w+') as resscript:
                                # only need first 11 lines of runscript.
                                with open(runscript_path, 'r') as runscript:
                                    scriptlines = runscript.readlines()[:12]
                                    scriptlines[7] = '#SBATCH -o {:s}\n'.format(target_of(outfile).replace('.o', 're.o'))
                                    scriptlines[8] = '#SBATCH -e {:s}\n'.format(target_of(outfile).replace('.o', 're.e'))
                                # In 'resume' script, insert command to resume the run using the last photo:
                                scriptlines.append('./re {:s}\n'.format(lastphoto))
                                resscript.writelines(scriptlines)

                            # Now make sure the desired stopping condition is placed in the inlist:
                            inlist_path = os.path.join(trackdir, 'inlist_project')
                            with open(inlist_path, 'r') as inlist:
                                inlistlines = inlist.readlines()

                            # add in new stopping condition and comment out the old one.
                            condition_infile = False
                            condition = 'stop_at_TP = .true.'
                            for j, inlistline in enumerate(inlistlines):
                                if '!gamma_center_limit' in inlistline:
                                    inlistlines[j] = inlistline.replace('!', '')
                               
                                if isuncommented_and_contains(inlistline, 'HB_limit', '!'):
                                    inlistlines[j] = inlistline.replace('HB_limit', '!HB_limit')
                                    # place the new stopping condition in after this line:
                                    insertion_idx = j+1
                                if 'stop_at_TP' in inlistline:
                                    condition_infile = True
                                    break

                            try:
                                if not condition_infile:
                                    # len(condition) + 7 has '7' b/c of 5 prepended spaces and two chars for '\n'
                                    inlistlines.insert(insertion_idx, ('{:s}\n'.format(condition)).rjust(len(condition)+7))
                            except NameError:
                                # 'insertion_index' was not defined.
                                print('{:s} did not have an uncommented \'HB_limit\' line.'.format(inlist_path))

                            # Now write the lines to a new inlist file:
                            newinlist_path = os.path.join(trackdir, 'inlist')
                            with open(newinlist_path, 'w+') as newinlist:
                                newinlist.writelines(inlistlines)
                            # replace old project_inlist too:
                            with open(inlist_path, 'r+') as inlist:
                                inlist.writelines(inlistlines)
                            # copy the re executable from MESAWORK_DIR.
                            shutil.copy(os.path.join(os.environ['MESAWORK_DIR'], 'cleanworkdir/re'), trackdir)
                            copy_tree(os.path.join(os.environ['MESAWORK_DIR'], 'cleanworkdir/src'), os.path.join(trackdir, 'src'))

        if len(badtracks) > 0:
            res_dict[rundir] = badtracks
            print("{:d}/{:d} tracks require extension in {:s}.".format(len(badtracks), len(trackdirs), target_of(rundir)))
        else:
            print("No tracks in {:s} require extension.".format(target_of(rundir)))
        print("===========================================================")

    return res_dict

def isuncommented_and_contains(string, phrase, comment_char='!'):

    comment_pos = string.find(comment_char)
    phrase_pos = string.find(phrase)
    
    # if (phrase is in srting) AND (e.g. the namelist comment '!' is not in the line OR comes after the phrase):
    return (phrase in string) & ( (comment_pos < 0) | (comment_pos > phrase_pos) )

def massof(trackname):

    # converts MIST track dir names into their resp. masses.

    mass_str = trackname.split('M')[0]

    return float(mass_str)/100.0

def target_of(path):
  
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
    mode = sys.argv[2]

    #prepend MIST grid path to runs:
    runname_str = os.path.join(os.environ['MIST_GRID_DIR'], 'MIST_v1.0', runname_str)

    if mode == 'rerun':
        # now check the runs for failed tracks:
        failed_dict = checkruns(runname_str)

        # Go through various run directories found...
        for rundir in failed_dict:
            # If there were failed tracks, re-submit them.
            failed_massdirs = failed_dict[rundir]
            print("___________________________________________________________")
            print("In {:s} re-submitting {:d} jobs.".format(target_of(rundir), len(failed_massdirs)))

            # Re-submit the failed tracks:
            for failed_massdir in failed_massdirs:
                os.chdir(failed_massdir)
                print("re-submitting {:s}...".format(target_of(failed_massdir)))
                slurmfile = glob(os.path.join(failed_massdir, '*M*run.sh'))[0]
                print(slurmfile)
                #os.system("sbatch " + slurmfile)
                os.chdir(os.environ['MIST_CODE_DIR'])

    elif mode == 'resume':
        res_dict = getphotos(runname_str)
        for rundir in res_dict:
            resume_dirs = res_dict[rundir]

            # Resume tracks found:
            print("___________________________________________________________")
            print("In {:s} resuming {:d} jobs.".format(target_of(rundir), len(resume_dirs)))

            # Re-submit the failed tracks:
            for resume_dir in resume_dirs:
                os.chdir(resume_dir)
                print("resuming {:s}...".format(target_of(resume_dir)))
                slurmfile = glob(os.path.join(resume_dir, '*M*resume.sh'))[0]
                #print(slurmfile)
                os.system("sbatch "+slurmfile)
                os.chdir(os.environ['MIST_CODE_DIR'])
                
