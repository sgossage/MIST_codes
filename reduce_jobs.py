#!/usr/bin/env python

"""

Takes MESA grids to create isochrones and organizes all the output files
into a nice directory structure.

The directory structure is as follows:
    top level directory --> MIST_vXX/FIDUCIAL/feh_pX.X_afe_pX.X/
    five subdirectories --> tracks/    eeps/    inlists/    isochrones/    plots/

Args:
    runname: the name of the grid

"""

import os
import sys
import subprocess

from scripts import mesa_plot_grid
from scripts import make_eeps_isos
from scripts import mesa_hist_trim
from scripts import reduce_jobs_utils
    
if __name__ == "__main__":
    
    runname = sys.argv[1]
    
    #Rename the run directory XXX as XXX_raw
    rawdirname = runname+"_raw"
    os.system("mv " + os.path.join(os.environ['MIST_GRID_DIR'],runname) + " " + os.path.join(os.environ['MIST_GRID_DIR'],rawdirname))
    
    #The XXX directory will contain the organized, reduced information
    newdirname = os.path.join(os.environ['MIST_GRID_DIR'],runname)
    os.mkdir(newdirname)
    
    #Make the eeps directory that will be filled in later
    os.mkdir(os.path.join(newdirname, "eeps"))
    
    #Make the isochrones directory that will be filled in later
    os.mkdir(os.path.join(newdirname, "isochrones"))

    print "************************************************************"
    print "****************SORTING THE HISTORY FILES*******************"
    print "************************************************************"
    reduce_jobs_utils.sort_histfiles(rawdirname)
    
    print "************************************************************"
    print "****************GENERATING A SUMMARY FILE*******************"
    print "************************************************************"
    reduce_jobs_utils.gen_summary(rawdirname)
    
    #Copy the summary file
    os.system("mv tracks_summary.txt " + newdirname)
    
    print "************************************************************"
    print "****************SORTING THE INLIST FILES********************"
    print "************************************************************"
    reduce_jobs_utils.save_inlists(rawdirname)
    
    print "************************************************************"
    print "***************SORTING THE PHOTOS AND MODELS****************"
    print "************************************************************"
    reduce_jobs_utils.save_lowM_photo_model(rawdirname)
    
    print "************************************************************"
    print "**********************MAKE ISOCHRONES***********************"
    print "************************************************************"
    #make_eeps_isos.make_eeps_isos(runname, basic=True)
    make_eeps_isos.make_eeps_isos(runname, basic=False)
    
    #print "************************************************************"
    #print "******************PLOTTING THE EEPS FILES*******************"
    #print "************************************************************"
    #os.mkdir(os.path.join(newdirname, "plots"))
    #mesa_plot_grid.plot_HRD(runname)
    #mesa_plot_grid.plot_combine(runname, iso=False, remove_pdf=False)

    #print "************************************************************"
    #print "******************PLOTTING THE ISOCHRONES*******************"
    #print "************************************************************"
    #mesa_plot_grid.plot_iso(runname)
    #mesa_plot_grid.plot_combine(runname, iso=True, remove_pdf=False)
    
    print "************************************************************"
    print "******COMPRESSING BOTH TRACKS AND REDUCED DIRECTORIES*******"
    print "************************************************************"
    #make a separate tracks directory
    os.system("mv " + os.path.join(newdirname, "tracks") + " " + newdirname + "_tracks")
    os.system("cp " + os.path.join(newdirname, "tracks_summary.txt") + " " + newdirname + "_tracks")

    os.chdir(os.environ['MIST_GRID_DIR'])    
    #When decompressed, this .tar.gz opens a MIST_vXX/feh_XXX_afe_XXX directory
    os.system("tar -zcvf " + '_'.join(runname.split('/')) + ".tar.gz " + runname)
    os.system("tar -zcvf " + '_'.join(runname.split('/')) + "_tracks.tar.gz " + runname+'_tracks')
    
    print "************************************************************"
    print "****************MIGRATING FILES TO STORAGE******************"
    print "************************************************************"
    os.system("rm -rf " + runname)
    os.system("rm -rf " + runname + '_tracks')
    print('Moving {:s} to {:s}...'.format(rawdirname, os.path.join(os.environ['STORE_DIR'], runname.split('/')[0])))
    os.system("mv " + rawdirname + " " + os.path.join(os.environ['STORE_DIR'], runname.split('/')[0]))
    print('Moving {:s} to {:s}...'.format('_'.join(runname.split('/')) + ".tar.gz ", os.path.join(os.environ['STORE_DIR'], runname.split('/')[0])))
    os.system("mv " + '_'.join(runname.split('/')) + ".tar.gz " + os.path.join(os.environ['STORE_DIR'], runname.split('/')[0]))
    print('Moving {:s} to {:s}...'.format('_'.join(runname.split('/')) + "_tracks.tar.gz ", os.path.join(os.environ['STORE_DIR'], runname.split('/')[0])))
    os.system("mv " + '_'.join(runname.split('/')) + "_tracks.tar.gz " + os.path.join(os.environ['STORE_DIR'], runname.split('/')[0]))

