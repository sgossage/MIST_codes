#!/usr/bin/env python

import os
import sys
import subprocess

from MIST_codes.scripts import make_eeps_isos

def remk_isos(runname):

    """
        This function will remake the .iso file from a grid run of MESA models using Jieun Choi's code (MIST) to create the .iso 
    files (which operates in conjuction with Aaron Dotter's iso code). This relies on the .tar.gz archives stored by the MIST codes
    after a grid run.

    runname (str):
    --------------
        This is the name of the run whose .iso file will be recreated. It should look something like:
  
        >> feh_mX.XX_afe_pX.X_vvcritX.X

        (maybe varied according to your naming convention when using MIST codes). MIST stores the tracks, eeps, and
    .iso files in a tar called <runname>.tar.gz in some specified storage directory by default.

    """

    # Need to untar the stored tracks(?)/eeps from a desired run.
    # Having done that, the .iso file can be recreated if told to
    # operate withing the directory that the files were untarred in.
    init_dir = os.getcwd()
    outdir_name = "test_output"

    storage_dir = '/n/conroyfs1/sgossage/MIST_grids/MIST_v1.0'#os.environ['STORE_DIR']
    output_dir = os.path.join(storage_dir, outdir_name)

    # Name corresp. to the tar containing the tracks of the run:
    #runname_tracks = runname + '_tracks'

    # Now need to untar the runname.tar.gz & the runname_tracks.tar.gz.
    # These will be untarred into the storage_dir/output directory. The
    # resulting structure will be:
    #
    #    storage_dir/output/MISTv1.0/runname
    #
    #    within lies: eeps/, isochrones/, inlists/, models_photos/
    #    once the runname_tracks.tar.gz is also untarred, this will
    #    contain tracks/ as well.
    #
    # The MIST code will be told to operate from this directory to recreate
    # the .iso file for the given run.

    # First task is to untar the run's .tar.gz archive so that it may be accessed:
    os.chdir(storage_dir)

    # The script below takes the runname as input and untars the contents within the output/
    # directory, creating the structure outlined above. It first untars the runname.tar.gz
    # and then the corresponding runname_tracks.tar.gz.

    print('Unpacking MIST_v1.0_{:s}.tar.gz (eeps, inlists, isos, photos)...'.format(runname))
    print(os.getcwd())
    os.system('./untar_run_to_output.sh ' + runname)
    #print('Unpacking MIST_v1.0_{:s}.tar.gz (the tracks)...'.format(runname_tracks))
    #os.system('./untar_run_to_output.sh ' + runname_tracks)
    print('...done!\n')

    # The tracks are presently stored in a directory called runname_tracks/, they will now
    # be moved into the runname/ directory:
    # One of my tar archives was structured in a strange way. The below segment checks the
    # structure and fixes it according to the peculiarity of my case, it is not a gneral fix.
    wrong_path = os.path.join(outdir_name, 'n')
    if os.path.isdir(wrong_path):
        right_path = os.path.join(wrong_path, 'regal', 'conroy_lab','sgossage','MIST_grids','*')
        os.system('mv {:s} {:s}'.format(right_path, outdir_name))
        os.system('rm -rf {:s}'.format(wrong_path))
    
    #track_destdir = os.path.join(outdir_name, 'MIST_v1.0', runname_tracks, '*')
    run_destdir =  os.path.join(outdir_name, 'MIST_v1.0', runname)

    # Make the directory that will house the run's tracks:
    #new_trackdir = os.path.join(run_destdir, 'tracks')
    #os.system('mkdir {:s}'.format(new_trackdir))

    # Move the tracks/remove the old (now empty directory from which they were moved):
    #print('Moving tracks into the same directory as the output...')
    #os.system('mv {:s} {:s}'.format(track_destdir, new_trackdir))
    #os.system('rm -rf {:s}'.format(track_destdir.split('*')[0]))
    #print('...done!\n')

    os.chdir(init_dir)

    # Now the files have been unpacked so we may run the MIST scripts in order to recreate
    # the .iso file for the run. The codes just need to be told to operate in the new location,
    # so this will involve modifying the make_eeps_isos script.

    MIST_runname = os.path.join('MIST_v1.0', runname)

    print('Creating the .iso file...\n')
    make_eeps_isos.make_eeps_isos(MIST_runname, basic=False, custom_path=output_dir)

    os.chdir(storage_dir)
    # Now transport the file to the output/ directory...
    # First make a backup of the .iso that is being replaced:
    print('\nBacking up old .iso file.')
    replace_dir = os.path.join('MIST_v1.0/output', runname, 'isochrones')
    iso_file_to_replace = os.path.join(replace_dir, 'MIST_v1.0_{:s}_full.iso'.format(runname))
    if os.path.isfile(iso_file_to_replace):
        os.system('cp {:s} {:s}.bckp'.format(iso_file_to_replace, iso_file_to_replace))
    else:
        print('Nevermind, it doesn\'t exist!')

    # Backup the old .cmd file:
    print('Backing up old .bin file.')
    if os.path.isfile("{:s}.bin".format(iso_file_to_replace)):
        os.system('mv {:s}.cmd {:s}.bin.bckp'.format(iso_file_to_replace, iso_file_to_replace))
    else:
        print('Nevermind, it doesn\'t exist!')

    # Now move the recreated .iso file to replace the old one:
    print('Moving new .iso file.')
    new_iso_file = os.path.join(outdir_name, 'MIST_v1.0', runname, 'isochrones', 'MIST_v1.0_{:s}_full.iso'.format(runname))
    os.system('mv {:s} {:s}'.format(new_iso_file, replace_dir))

    # Finally, remove the unpacked files, they are unnecessary.
    print('Deleting untarred files (.tar.gz still intact).\n')
    os.system('rm -rf {:s}'.format(os.path.join(output_dir, '*')))
    os.chdir(init_dir)

    return

if __name__ == '__main__':

    runname = sys.argv[1]

    remk_isos(runname)
