import os
import glob

def createcmd(iso_filename, photstr='UBVRIplus', Av=None, gravdark_i = 0.0, force=False):

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

    filename = iso_filename

    # Acquire the .iso files:
    # Directory to A. Dotter's isochrone code on my system:
    isocodedir = os.environ['ISO_DIR']

    # Switch to the iso code directory:
    os.chdir(isocodedir)

    out_fname = filename.split('.iso')[0] +'_{:s}'.format(photstr) + '.iso.cmd'
    # check if the .cmd file already exists; it may not need to be recreated:
    if os.path.isfile(out_fname) and not force:
        #print("The .cmd file \"{:s}\" already exists...".format(out_fname))
        with open(out_fname, 'r') as cmdf:
            cmdflines = cmdf.readlines()
        # Check if extinction in the currently existing .cmd file matches desired value:
        if "{:.3f}".format(float(cmdflines[8].split()[-1])) == "{:.3f}".format(Av):
            #print("Desired Av = {:.3f} matches existing Av = {:.3f}; the existing .cmd file will be used.".format(Av, float(cmdflines[8].split()[-1])))
            pass
        else:
            #print("Desired Av = {:.3f} does not match existing Av = {:.3f}; creating a new .cmd file...".format(Av, float(cmdflines[8].split()[-1])))
            os.system('./make_cmd ' + photstr + ' ' + filename + ' ' + str(Av))

    else:
        print("Creating .cmd file: {:s}".format(out_fname))
        # For each .iso file, run the iso code:
        if Av == None:
            os.system('./make_cmd ' + photstr + ' ' + filename)
        else:
            os.system('./make_cmd ' + photstr + ' ' + filename + ' ' + str(Av))

    # Return to former directory:
    os.chdir(initdir)

    # output .cmd file is the same as the .iso file, but with .cmd appended to it:
    #out_fnames = [isofile + '.cmd' for isofile in isofiles]
    return out_fname
