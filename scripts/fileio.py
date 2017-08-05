import glob
import os
import numpy as np
import matplotlib.pyplot as plt

import isointerp
import createcmd as ccmd
import isomist
from gdiso import *

def get_fn(feh, vvcrit, mode, mass=0.0, ebv = 0.0, gravdark_i = 0.0, exttag=None):
    """
        Gets the filename, given a [Fe/H] value & a v/vcrit value. Assumes that the environment variable 'STORE_DIR' is set
        and that the grid names follow feh_pX.XX_afe_0.0_vvcritX.X for their name format.

        Args:
            feh: the [Fe/H] of the desired grid from which the cmd filename will be taken.
            vvcrit: same as above, but for v/vcrit
            mode: a string -- either 'iso' for .iso files or 'eep' for .track.eep files.
            ebv: E(B-V) desired for the .cmd file.
            gravdark_i: inclination angle for gravity darkening; default is 0.0 angles, i.e. looking @ the equator of the stars.
            exttag: any tags to append to the end of the grid; e.g. if feh_pX.XX_afe_0.0_vvcritX.X_HB was desired, this would be 'HB'
                    (a str, & '_' is added by this code itself).
    """

    store_dir =  os.path.join(os.environ['STORE_DIR'], 'MIST_v1.0', 'output')
    if feh < 0.0:
        fehstr = 'm{:.2f}'.format(abs(feh))
    else:
        fehstr = 'p{:.2f}'.format(feh)

    rotstr = '{:.1f}'.format(vvcrit)

    # Getting the track's file name:
    if exttag == None:
        gridname = 'feh_{:s}_afe_p0.0_vvcrit{:s}'.format(fehstr, rotstr)
    else:
        gridname = 'feh_{:s}_afe_p0.0_vvcrit{:s}_{:s}'.format(fehstr, rotstr, exttag)

    grid_dir = os.path.join(store_dir, gridname)

    if mode == 'iso' or mode == 'isocmd':
        # Av extinction value
        Av = 3.1*ebv

        if mode == 'iso':
            isof = get_isocmdf(grid_dir, feh, vvcrit, Av = Av, gravdark_i = gravdark_i, exttag = exttag, create_cmd = False)
            #assert len(isof) == 1, "Got {:d}, rather than one .cmd file in {:s}.".format(len(isof), grid_dir)
            return isof
        elif mode == 'isocmd':
            # get the .iso.cmd file name:
            cmdf = get_isocmdf(grid_dir, feh, vvcrit, Av = Av, gravdark_i = gravdark_i, exttag = exttag, create_cmd = True)
            # in case more than one .cmd file is returned -- should never happen since filenames are unique; indicates
            # something may be wrong with file seaerch terms.
            #assert len(cmdf) == 1, "Got {:d}, rather than one .cmd file in {:s}.".format(len(cmdf), grid_dir)
            return cmdf


    if mode == 'eep':

        eepf = get_masstrackeepf(grid_dir, mass)
        # Trying to get the track.eep for a given mass value. Print masses found on failure:
        if not eepf:
            # if a track.eep file was not able to be retrieved (due to unavailable mass?), list the available masses here:
            filelist = glob.glob(os.path.join(grid_dir, 'eeps/*M.track.eep'))
            print('Valid masses found:\n')
            masslist = []
            for f in filelist:
                starmass = float(f.split('eeps/')[1].split('M')[0])/100.0
                masslist.append(starmass)
            return sorted(masslist)
        # or if successful the recovered .track.eep file name will be returned.
        return eepf

def get_isocmdf(grid_dir, feh, vvcrit, Av = 0.0, gravdark_i=0.0, exttag=None, create_cmd=True):

    """
        Gets the .iso.cmd file for the fiven feh and v/vcrit values. If the file DNE, tjos wo;; attempt to 
    call A. Dotter's iso interpolation code in order to create it from existing .iso files.
    """

    # I don't think that filelist needs to be a list...it is meant to just represent one file.

    if os.path.isdir(grid_dir):

        # The .cmd file may still need to be created.
        if exttag == None:
            # checks for basic.iso first:
            isofile = glob.glob(os.path.join(grid_dir, 'isochrones/*basic.iso'))
            if len(isofile) == 0:
                isofile = glob.glob(os.path.join(grid_dir, 'isochrones/*full.iso'))

            assert len(isofile) == 1, "Got {:d}, rather than one .iso file in {:s}.".format(len(isofile), grid_dir)
            filename = isofile[0]
        else:
            isofile = glob.glob(os.path.join(grid_dir, 'isochrones/*{:s}_basic.iso'.format(exttag)))
            if len(isofile) == 0:
                isofile = glob.glob(os.path.join(grid_dir, 'isochrones/*full.iso'))

            assert len(isofile) == 1, "Got {:d}, rather than one .iso file in {:s}.".format(len(isofile), grid_dir)
            filename = isofile[0]

        print('-----')
        if create_cmd:
            cmdfile = isomist.createcmd(filename, Av = Av, gravdark_i = gravdark_i)
            filename = cmdfile

        return filename

    else:
        # If non-existent in the default grid, attempt interpolation.
        interpiso_fn = isointerp.call_isointerp(vvcrit, feh)
        filename = interpiso_fn
        if create_cmd:
            cmdfile = isomist.createcmd(interpiso_fn, Av = Av, gravdark_i = gravdark_i)
            filename = cmdfile

        return filename

def get_masstrackeepf(grid_dir, mass):

    """
        Gets the track.eep file for the given mass.
    """

    mass_string = ''.join(str(mass/100.0).split('.'))
    # The mass string is 5 digits...e.g. 8.0 = 00800, 0.3 = 00030. Add a zero for missing digits:
    if len(mass_string) < 5:
        for i in range(5 - len(mass_string)):
            mass_string = mass_string + '0'

    filelist = glob.glob(os.path.join(grid_dir, 'eeps/*{:s}M.track.eep'.format(mass_string)))
    if not filelist:
        print('No track.eep file was found for {:.1f} Msun.'.format(mass))
        return filelist
    else:
        return filelist[0]
