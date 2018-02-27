import glob
import os
import numpy as np
from MIST_codes.scripts.gdiso import *

def createcmd(isoobj, photstr='UBVRIplus', feh= None, vvcrit = None, cov = None, Av=None, fname = None, gravdark_i = 0.0):

    """
        This function will find a number of MIST .iso files and run A. Dotter's make_cmd code using them to create
    corresponding .cmd files, containing magnitudes, etc.

    From A. Dotter's iso code:
    ============================================================================================================
    make_cmd:   
    usage: ./make_cmd [phot string] [isochrone file] [Av]
      [phot string] = UBVRIJHKs, etc.                    
      [isochrone file] = name of isochrone file to transform
      [Av] optional argument; if not set then take value from input.nml
      all other options set through cmd_controls in input.nml
    ============================================================================================================
        So a photometric strng, the .iso file, and potentially an Av value should be supplied. Av enters in this
    function (createcmd()) as an optional input parameter, while the photometric string is a required input param.

    """
    initdir = os.getcwd()

    # Acquire the .iso files:

    # if a feh or vvcrit value is given, form a search string:
    search_str = ''
    if feh != None or vvcrit != None:
        if feh >= 0.00:
            search_str = search_str + '*_feh_p{:.2f}_'.format(feh)
        elif feh < 0.00:
            search_str = search_str + '*_feh_m{:.2f}_'.format(abs(feh))
        if cov is not None and cov != 0.016:
            search_str = search_str + "*vvcrit{:.1f}_cov{:.3f}_full".format(vvcrit, cov)
        else:
            search_str = search_str + '*vvcrit{:.1f}_full'.format(vvcrit)

        search_str = search_str + '.iso'
    #print("!!!", search_str) 

    # Directory leading to the stored .iso files (for MY environment):
    storedir = os.environ['STORE_DIR']
    if fname == None:
        pathtoisos = os.path.join(storedir, 'MIST_v1.0/output/*/isochrones/{:s}'.format(search_str))
    else:
        pathtoisos = fname
    #print(pathtoisos)
    # Returns a list of paths to each .iso file found in path structure above:
    isofiles = glob.glob(pathtoisos)
    # Go through the .iso file list and run A. Dotter's make_cmd program on them:
    # Directory to A. Dotter's isochrone code on my system:
    isocodedir = os.environ['ISO_DIR']

    # Switch to the iso code directory:
    os.chdir(isocodedir)

    out_fnames = []
    for idx, filename in enumerate(isofiles):

        # Check if gravity darkening is desired:
        if gravdark_i > 0.0:
            # Will check for GDed .iso file; creates one if nec.
            isoobj = ISO(filename)
            filename = rw_gdiso(isoobj, filename, gravdark_i)
            # replace filename with the new GD file:
            isofiles[idx] = filename

        out_fnames.append(filename + '.cmd')
        # check if the .cmd file already exists; it may not need to be recreated:
        if os.path.isfile(out_fnames[idx]):
            #print("The .cmd file \"{:s}\" already exists...".format(out_fnames[idx]))
            with open(out_fnames[idx], 'r') as cmdf:
                cmdflines = cmdf.readlines()
            # Check if extinction in the currently existing .cmd file matches desired value:
            if "{:.3f}".format(float(cmdflines[8].split()[-1])) == "{:.3f}".format(Av):
                # if it does, continue to next file
                #print("Desired Av = {:.3f} matches existing Av = {:.3f}; the existing .cmd file will be used.".format(Av, float(cmdflines[8].split()[-1])))
                continue
            else:
                #print("Desired Av = {:.3f} does not match existing Av = {:.3f}; creating a new .cmd file...".format(Av, float(cmdflines[8].split()[-1])))
                # or else, continue this iteration and create a new .cmd file:
                pass

        # For each .iso file, run the iso code:
        #print filename
        if Av == None:
            #agemax = np.amax(np.genfromtxt(filename, usecols=(1,), skip_header = 11))
            #agemin = np.amin(np.genfromtxt(filename, usecols=(1,), skip_header = 11))
            #print('\nCreating .cmd file from {:s} (log10 age range {:f} to {:f})...\n'.format(filename,agemin,agemax))
            os.system('./make_cmd ' + photstr + ' ' + filename)
        else:
            #agemax = np.amax(np.genfromtxt(filename, usecols=(1,), skip_header = 11))
            #agemin = np.amin(np.genfromtxt(filename, usecols=(1,), skip_header = 11))
            #print('\nCreating .cmd file from {:s} (log10 age range {:f} to {:f})...\n'.format(filename,agemin, agemax))
            os.system('./make_cmd ' + photstr + ' ' + filename + ' ' + str(Av) + ' ' + '.false.')

    # Return to former directory:
    os.chdir(initdir)

    # output .cmd file is the same as the .iso file, but with .cmd appended to it:
    #out_fnames = [isofile + '.cmd' for isofile in isofiles]
    return out_fnames

    # Reporting the max/min ages of the .cmd file passed back (debugging):
    #for output_pathname in out_fnames:
    #    agemax = np.amax(np.genfromtxt(output_pathname, usecols=(1,), skip_header = 13))
    #    agemin = np.amin(np.genfromtxt(output_pathname, usecols=(1,), skip_header = 13))
    #    print('\nThe .cmd file, {:s}, has age range: {:f} to {:f} in log10 age.'.format(output_pathname, agemin, agemax))
    #    print('\n')

    # Now we are done, the .cmd files should have been created in the same directories as the .iso files.

    #if gravdark_on:
        # Switch on conversion of L/Teff to mimick grav. darkening:
    #    print('\nEnabling gravity darkening effect.\n')
    #    with open('input.nml', 'r') as infile:
    #        lines = infile.readlines()
    #        for j, line in enumerate(lines):
    #            if 'include_gravity_darkening' in line:
    #                lines[j] = line.replace('.false.', '.true.')
    #    with open('input.nml', 'w') as outfile:
    #        for line in lines:
    #            outfile.write(line)

    #elif not gravdark_on:
    #    print('\nDisabling gravity darkening effect.\n')
    #    with open('input.nml', 'r') as infile:
    #        lines = infile.readlines()
    #        for j, line in enumerate(lines):
    #            if 'include_gravity_darkening' in line:
    #                lines[j] = line.replace('.true.', '.false.')
    #    with open('input.nml', 'w') as outfile:
    #        for line in lines:
    #            outfile.write(line)
