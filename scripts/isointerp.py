import glob
import os
import numpy as np
from pathing import isopath

def call_isointerp(vvcrit, feh, cov=None, custom_isopath = None, exttag = None):

    """
        In order to get closer to a calcsfh best fit, may need to interpolate available MIST isochrones. A. Dotter's
    iso code can do this, if supplied with a list of available isochrone files. The vvcrit value will single out the group
    of metallicities to use when interpolating.

    Order matters in the logZ list: go from lower to higher values (e.g., -0.30, ..., 0.30); this function handles that.

    This returns the path to the interpolated .iso file.

    """
    # Exit if the metallicity should already exist (not a great security algorithm atm?):
    feh_masterlist = [-0.75, -0.60, -0.45, -0.30, -0.15, 0.00, 0.15, 0.30, 0.45]
    if feh in feh_masterlist and custom_isopath == None:
        print('Metallicity value given already exists; no need to interpolate.')
        return

    # Track the initial directory:
    initdir = os.getcwd()

    # Get the available isochrone files. These are stored in /n/conroyfs1/.../output/:
    if custom_isopath == None:
        isofiles = isopath(vvcrit, cov, exttag=exttag)
    else:
        isofiles = glob.glob(os.path.join(custom_isopath, '*.iso'))

    isodict_list = []
    fehdict_list = []
    for afile in isofiles:
        fehstr, fehval = fname_getfeh(afile)
        isodict_list.append({afile: fehval})
        fehdict_list.append({fehstr: fehval})

    # We now have two lists of dictionaries mapping the file names to their feh values, and the feh values to feh strings.
    # These may be sorted in ascending order like so...
    # First, compile the list of dictionaries into a single dictionary:
    isodict = compile_dictlist(isodict_list)
    fehdict = compile_dictlist(fehdict_list)

    # We will sort the dictionary keys (file names or feh strings) by their values (feh values):
    isofname_list = sorted(isodict, key=isodict.__getitem__)
    #print(isofname_list)
    fehstr_list = sorted(fehdict, key=fehdict.__getitem__)
    #print(fehstr_list)
    filestr_lines = [(fehstr_list[i], isofname_list[i]) for i in range(len(feh_masterlist))]

    # isofiles is now a list of the .iso files available at the given vvcrit value.
    # The iso code requires a list of these files to interpolate isochrones, so let's make it:

    # The format is:
    #----------------------------------
    #  <Num iso files>
    #  <file name>
    #  ...
    #  <file name>
    #----------------------------------
    # The list should go in ascending order of feh values from top to bottom. 

    isocodedir = os.environ['ISO_DIR'] 
    # Open the list file for writing:
    interplist_fname = os.path.join(isocodedir, 'interp.list')
    interplist_f = open(interplist_fname, 'w')
    # Now write according to the format:
    interplist_f.write("{:d}\n".format(len(feh_masterlist)))
    for i in range(len(filestr_lines)):
        interplist_f.write("{:s}\n".format(filestr_lines[i][1]))
    
    # and that should be all that is necessary to create the interpolation file list for the iso code. The list is in the iso code directory
    # under the name 'interp.list'.
    interplist_f.close()

    # Now call A. Dotter's code to perform the interpolation:
    #  iso_interp_met        
    #    usage:             
    #    ./iso_interp_met [list] [Z/H] [output]

    os.chdir(isocodedir)

    # Make the output file name; using the path to interpolated files relevant to my system:
    interpolation_dir = os.path.join(os.environ['STORE_DIR'], 'MIST_v1.0', 'output', 'interpolations')#'/home/seth/Research/MIST/MISTout/MIST_v1.0/output/interpolations'

    if feh >= 0.00:
        if exttag is not None:
            output_fname = 'MIST_v1.0_feh_p{:.2f}_afe_p0.0_vvcrit{:.1f}_{:s}_interp.iso'.format(feh, vvcrit, exttag)
        else:
            output_fname = 'MIST_v1.0_feh_p{:.2f}_afe_p0.0_vvcrit{:.1f}_interp.iso'.format(feh, vvcrit)
    elif feh < 0.00:
        if exttag is not None:
            output_fname = 'MIST_v1.0_feh_m{:.2f}_afe_p0.0_vvcrit{:.1f}_{:s}_interp.iso'.format(abs(feh), vvcrit, exttag)
        else:
            output_fname = 'MIST_v1.0_feh_m{:.2f}_afe_p0.0_vvcrit{:.1f}_interp.iso'.format(abs(feh), vvcrit)

    output_pathname = os.path.join(interpolation_dir, output_fname)    

    # Call the code:
    if os.path.isfile(output_pathname):
        pass
    else:
        print('Creating .iso file: {:s}.'.format(output_pathname))
        os.system('./iso_interp_met interp.list {:f} {:s}'.format(feh, output_pathname))

    os.chdir(initdir)
    
    # Now the interpolated .iso file should have been created.
    #print(output_pathname)
    agemax = np.amax(np.genfromtxt(output_pathname, usecols=(1,), skip_header = 11))
    agemin = np.amin(np.genfromtxt(output_pathname, usecols=(1,), skip_header = 11))
    print('\nThe created .iso file, {:s}, has age range: {:f} to {:f} in log10 age.\n'.format(output_pathname, agemin, agemax))

    return output_pathname

def fname_getfeh(fname):

    """
        This returns the feh value in a file name according to my naming scheme: *_feh_<fehstring>_afe_<afestring>_vvcrit_<vvcritstring>_*.
    """

    fname = os.path.basename(fname)

    fehstr = (fname.split('feh_')[-1]).split('_afe')[0]

    if 'p' in fehstr:
        if float(fehstr.split('p')[-1]) == 0.00:
            fehstr = fehstr.split('p')[-1]
        elif float(fehstr.split('p')[-1]) > 0.00: 
            fehstr = '+' + fehstr.split('p')[-1]

    elif 'm' in fehstr:
        fehstr = '-' + fehstr.split('m')[-1]

    # The actual feh value as a floating point number:
    feh = float(fehstr)    

    return fehstr, feh

def compile_dictlist(dictlist):

    """
        This goes through a list of dictionaries and compiles each member into a single dictionary.
    """

    compdict = { key: value for dictionary in dictlist for key, value in dictionary.items() }

    # Return the compiled dictionary:
    return compdict
