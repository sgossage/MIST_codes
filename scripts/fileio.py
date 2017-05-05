import glob
import os
import numpy as np
import matplotlib.pyplot as plt

import isointerp
import createcmd as ccmd
from gdiso import *


def get_isocmdf(grid_dir, feh, vvcrit, Av = 0.0, gravdark_i=0.0, exttag=None):

    """
        Gets the .iso.cmd file for the fiven feh and v/vcrit values. If the file DNE, tjos wo;; attempt to 
    call A. Dotter's iso interpolation code in order to create it from existing .iso files.
    """

    # I don't think that filelist needs to be a list...it is meant to just represent one file.

    if os.path.isdir(grid_dir):

        # The .cmd file may still need to be created.
        if exttag == None:
            isofile = glob.glob(os.path.join(grid_dir, 'isochrones/*full.iso'))[0]
        else:
        	isofile = glob.glob(os.path.join(grid_dir, 'isochrones/*{:s}_full.iso'.format(exttag)))[0]
    
        # Call Aaron's orientation code to recompute Teff and L & create a GDed .cmd
        filelist = ccmd.createcmd(fname = isofile, Av = Av, gravdark_i = gravdark_i)
                
        return filelist

    else:
        # If non-existent in the default grid, attempt interpolation.
        interp_isof = isointerp.call_isointerp(vvcrit, feh)
        filelist = ccmd.createcmd(fname = interp_isof, Av = Av, gravdark_i = gravdark_i)
        return filelist


def plot_photof(fname, erron= False, ax=None):

    """
        Plots the magnitudes stored in a MATCH .phot file. This is no a special file, aside from the format of no header,
    and only two columns. Column 1 = bluer magniudes, column 2 = redder magnitudes; this function extracts the magnitudes
    and plots them on a CMD of red vs. blue - red.
    """

    if not ax:
        ax = plt.gca()

    V = np.genfromtxt(fname, usecols = (0, ))
    I = np.genfromtxt(fname, usecols = (1, ))
    colors = V - I
    x = colors
    y = V

    ax.scatter(x, y, c='b', label="{:s}".format(fname.split('/')[-1].split('.phot')[0]), s = 2,  alpha = 0.6, zorder = 9999)

    if erron:
        errfile = fname.split('.phot')[0] + '.err'

        if os.path.isfile(errfile):

            verr = np.genfromtxt(errfile, usecols = (0, ))
            ierr = np.genfromtxt(errfile, usecols = (1, ))
            color_err = np.sqrt(verr**2 + ierr**2)
            xerr = color_err
            yerr = verr

            ax.errorbar(x, y , xerr=xerr, yerr=yerr, fmt='o', marker=None, alpha = 0.4)

    return ax

def plot_geneva_iso(lage, vvcrit, pltmass, plot_index, color_n, ax):

    """

        Not very flexible code at the moment. Only works for Geneva models that I have downloaded and
    these only exist at solar Z, lage = 8.2, and v/vcrit = 0.0, 0.6.

    This expects a naming convention of:
    
        Isochrones_lage<val>_vvcrit<val>_solarZ.dat

    for the Geneva isochrone file. This function extracts the log L and log Teff information and plots it.
    I have these files store in my home directory in a folder called "Geneva" and these files are sought
    from that location.

    """

    geneva_dir = '/n/conroyfs1/sgossage/Geneva'

    geneva_isos = glob.glob(os.path.join(geneva_dir, '*.dat'))
    #print(geneva_isos)
    for isoname in geneva_isos:
        isoname = isoname.split('/')[-1]
        if "lage{:.2f}".format(lage) in isoname and "vvcrit{:.1f}".format(vvcrit) in isoname:
            # If the isochrone exists, plot it. Right now only does log L and log Teff,
            # but maybe also do filters...some are available in the Geneva isochrone file,
            # although it will need to match the MIST isochrone file's filter set.
            logL = np.genfromtxt(os.path.join(geneva_dir, isoname),
                          usecols = (4,), skip_header = 2)

            logTeff = np.genfromtxt(os.path.join(geneva_dir, isoname),
                            usecols = (5,), skip_header = 2)

            init_masses = np.genfromtxt(os.path.join(geneva_dir, isoname),
                            usecols = (0,), skip_header = 2)

            x = logTeff
            y = logL

            linecolor = plt.cm.rainbow(color_n[plot_index])

            base_line, = ax.plot(x, y, ls = '--', c = linecolor, label='Geneva: age = {:.2e}, [Fe/H] = 0.00, vvcrit = {:.1f}'.format(10**lage, vvcrit), 
                                 alpha = 0.8, zorder=-9999)

            print("Plotting {:s}".format(isoname))

            # pltmass may be given as a list of floats in order to plot a specific mass.
            if isinstance(pltmass, float):
                # Need to get the appropriae index of the given mass....should search for closest mass.
                diff_arr = abs(init_masses - pltmass)
                m_i = np.where(diff_arr == np.min(diff_arr))[0][0]
                print('Plotting mass {:.2f} Msun ({:.2f} Msun requested).'.format(init_masses[m_i], pltmass))
                ax.scatter(x[m_i], y[m_i], lw=0.1, alpha=0.5, color = base_line.get_color(), zorder=2)
                ax.text(x[m_i]*0.98, y[m_i], '{:.2f}'.format(init_masses[m_i]) + r' $M_{\odot}$', fontsize=4, color = base_line.get_color())

            elif isinstance(pltmass, list):
                for mass in pltmass:
                    # Need to get the appropriae index of the given mass....should search for closest mass.
                    diff_arr = abs(init_masses - mass)
                    m_i = np.where(diff_arr == np.min(diff_arr))[0][0]
                    print('Plotting mass {:.2f} Msun ({:.2f} Msun requested).'.format(init_masses[m_i], mass))
                    ax.scatter(x[m_i], y[m_i], lw=0.1, alpha=0.5, color = base_line.get_color(), zorder=2)
                    ax.text(x[m_i]*0.98, y[m_i], '{:.2f}'.format(init_masses[m_i]) + r' $M_{\odot}$', fontsize=4, color = base_line.get_color())

            break 

    return