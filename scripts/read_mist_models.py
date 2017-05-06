from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib as mpl
mpl.use('Agg')
import os
import glob
import seaborn

from MIST_codes.scripts.fileio import *

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ISO:
    
    """
    
    Reads in MIST isochrone files.

    
    """
    
    def __init__(self, filename, verbose=True):
    
        """
        
        Args:
            filename: the name of .iso file.
        
        Usage:
            >> iso = read_mist_models.ISO('MIST_v1.0_feh_p0.00_afe_p0.0_vvcrit0.4.iso')
            >> age_ind = iso.age_index(8.0)
            >> logTeff = iso.isos[age_ind]['log_Teff']
            >> logL = iso.isos[age_ind]['log_L']
            >> plt.plot(logTeff, logL) #plot the HR diagram for logage = 8.0
            
        Attributes:
            version     Dictionary containing the MIST and MESA version numbers.
            abun        Dictionary containing Yinit, Zinit, [Fe/H], and [a/Fe] values.
            rot         Rotation in units of surface v/v_crit.
            ages        List of ages.
            num_ages    Number of isochrones.
            hdr_list    List of column headers.
            isos        Data.
            
        """
        
        self.filename = filename
        if verbose:
            print('Reading in: ' + self.filename)
            
        self.version, self.abun, self.rot, self.ages, self.num_ages, self.hdr_list, self.isos = self.read_iso_file()
        
    def read_iso_file(self):

        """

        Reads in the isochrone file.
        
        Args:
            filename: the name of .iso file.
        
        """
        
        #open file and read it in
        with open(self.filename) as f:
            content = [line.split() for line in f]
        version = {'MIST': content[0][-1], 'MESA': content[1][-1]}
        abun = {content[3][i]:float(content[4][i]) for i in range(1,5)}
        rot = float(content[4][-1])
        num_ages = int(content[6][-1])
        
        #read one block for each isochrone
        iso_set = []
        ages = []
        counter = 0
        data = content[8:]
        for i_age in range(num_ages):
            #grab info for each isochrone
            num_eeps = int(data[counter][-2])
            num_cols = int(data[counter][-1])
            hdr_list = data[counter+2][1:]
            formats = tuple([np.int32]+[np.float64 for i in range(num_cols-1)])
            iso = np.zeros((num_eeps),{'names':tuple(hdr_list),'formats':tuple(formats)})
            #read through EEPs for each isochrone
            for eep in range(num_eeps):
                iso_chunk = data[3+counter+eep]
                iso[eep]=tuple(iso_chunk)
            iso_set.append(iso)
            ages.append(iso[0][1])
            counter+= 3+num_eeps+2
        return version, abun, rot, ages, num_ages, hdr_list, iso_set  
        
    def age_index(self, age):
    
        """

        Returns the index for the user-specified age.
    
        Args:
            age: the age of the isochrone.
    
        """
    
        diff_arr = abs(np.array(self.ages) - age)
        age_index = np.where(diff_arr == min(diff_arr))[0][0]
    
        if ((age > max(self.ages)) | (age < min(self.ages))):
            print('The requested age is outside the range. Try between ' + str(min(self.ages)) + ' and ' + str(max(self.ages)))
        
        return age_index

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    		
class ISOCMD:
    
    """
    
    Reads in MIST CMD files.

    
    """
    
    def __init__(self, feh, vvcrit, ebv=0.0, gravdark_i=0.0, exttag=None, filename=None, verbose=True):
    
        """
        
        Args:
            filename: the name of .iso.cmd file.
        
        Usage:
            >> isocmd = read_mist_models.ISOCMD('MIST_v1.0_feh_p0.00_afe_p0.0_vvcrit0.4.iso.cmd')
            >> age_ind = isocmd.age_index(7.0)
            >> B = isocmd.isocmds[age_ind]['Bessell_B']
            >> V = isocmd.isocmds[age_ind]['Bessell_V']
            >> plt.plot(B-V, V) #plot the CMD for logage = 7.0
        
        Attributes:
            version         Dictionary containing the MIST and MESA version numbers.
            photo_sys       Photometric system. 
            abun            Dictionary containing Yinit, Zinit, [Fe/H], and [a/Fe] values.
            Av_extinction   Av for CCM89 extinction.
            rot             Rotation in units of surface v/v_crit.
            ages            List of ages.
            num_ages        Number of ages.
            hdr_list        List of column headers.
            isocmds         Data.
        
        """

        self.filename = self.get_fn(feh, vvcrit, ebv = ebv, gravdark_i = gravdark_i, exttag=exttag)
        self.exttag = exttag
        self.feh = feh
        self.gdark_i = gravdark_i
        self.xextent = [0,0]
        self.yextent = [0,0]
 
        if verbose:
            print('Reading in: ' + self.filename)
            
        self.version, self.photo_sys, self.abun, self.Av_extinction, self.rot, self.ages, self.num_ages, self.hdr_list, self.isocmds = self.read_isocmd_file()
    
    def read_isocmd_file(self):

        """

        Reads in the cmd file.
        
        
        Args:
            filename: the name of .iso.cmd file.
        
        """
        
        #open file and read it in
        with open(self.filename) as f:
            content = [line.split() for line in f]
        version = {'MIST': content[0][-1], 'MESA': content[1][-1]}
        photo_sys = ' '.join(content[2][4:])
        abun = {content[4][i]:float(content[5][i]) for i in range(1,5)}
        rot = float(content[5][-1])
        num_ages = int(content[7][-1])
        Av_extinction = float(content[8][-1])
        
        #read one block for each isochrone
        isocmd_set = []
        ages = []
        counter = 0
        data = content[10:]
        for i_age in range(num_ages):
            #grab info for each isochrone
            num_eeps = int(data[counter][-2])
            num_cols = int(data[counter][-1])
            hdr_list = data[counter+2][1:]
            formats = tuple([np.int32]+[np.float64 for i in range(num_cols-1)])
            isocmd = np.zeros((num_eeps),{'names':tuple(hdr_list),'formats':tuple(formats)})
            #read through EEPs for each isochrone
            for eep in range(num_eeps):
                isocmd_chunk = data[3+counter+eep]
                isocmd[eep]=tuple(isocmd_chunk)
            isocmd_set.append(isocmd)
            ages.append(isocmd[0][1])
            counter+= 3+num_eeps+2
        return version, photo_sys, abun, Av_extinction, rot, ages, num_ages, hdr_list, isocmd_set

    def age_index(self, age):
        
        """

        Returns the index for the user-specified age.
        
        Args:
            age: the age of the isochrone.
        
        """
        
        diff_arr = abs(np.array(self.ages) - age)
        age_index = np.where(diff_arr == min(diff_arr))[0][0]
        
        if ((age > max(self.ages)) | (age < min(self.ages))):
            print('The requested age is outside the range. Try between ' + str(min(self.ages)) + ' and ' + str(max(self.ages)))
            
        return age_index

    def get_fn(self, feh, vvcrit, ebv = 0.0, gravdark_i = 0.0, exttag=None):
        """
            Gets the filename, given a [Fe/H] value & a v/vcrit value. Assumes that the environment variable 'STORE_DIR' is set
            and that the grid names follow feh_pX.XX_afe_0.0_vvcritX.X for their name format.

            Args:
                feh: the [Fe/H] of the desired grid from which the cmd filename will be taken.
                vvcrit: same as above, but for v/vcrit
                ebv: E(B-V) desired for the .cmd file.
                gravdark_i: inclination angle for gravity darkening; default is 0.0 angles, i.e. looking @ the equator of the stars.
                exttag: any tags to append to the end of the grid; e.g. if feh_pX.XX_afe_0.0_vvcritX.X_HB was desired, this would be 'HB'
                        (a str, & '_' is added by this code itself).
        """

        store_dir =  os.path.join(os.environ['STORE_DIR'])
        print(store_dir)
        if feh < 0.0:
            fehstr = 'm{:.2f}'.format(abs(feh))
        else:
            fehstr = 'p{:.2f}'.format(feh)

        rotstr = '{:.1f}'.format(vvcrit)

        # Get the .iso.cmd file name/path
        gridname = 'feh_{:s}_afe_p0.0_vvcrit{:s}'.format(fehstr, rotstr)
        print(gridname)
        # if there is an extra tag to append (e.g. for diff conv. core ovsh values):
        if exttag != None:
            gridname += '_{:s}'.format(exttag)

        grid_dir = os.path.join(store_dir, gridname)

        print(grid_dir)
        Av = 3.1*ebv

        cmdf = get_isocmdf(grid_dir, feh, vvcrit, Av = Av, gravdark_i = gravdark_i, exttag = exttag)
        assert len(cmdf) == 1, "Got more than one .cmd file in {:s}.".format(grid_dir)
        #self.filename = cmdf[0]
        return cmdf[0]

    def set_isodata(self, lage, x_name, y_name, dmod=0.0, ax=None, phasemask=[], x_to_ymx=True, geneva_on=False):
        """
            Given a log10 age, get the data desired to plot the CMD or HRD
        """
        exttag = self.exttag

        # For labeling age in potential plots:
        if lage >= 9.0:
            lage_str = "{:.1f} Gyr".format((10**lage)/(10.**9))
        if 6.0 <= lage < 9.0:
            lage_str = "{:.1f} Myr".format((10**lage)/(10.**6))

        if exttag != None:
            self.lbl = 'MIST({:s}): age = {:s} Myr, [Fe/H] = {:.2f}, vvcrit = {:.1f} (i = {:.1f})'.format(exttag, lage_str,  self.feh, self.rot, self.gdark_i*(180./np.pi))
        else:
            self.lbl = 'MIST: age = {:s}, [Fe/H] = {:.2f}, vvcrit = {:.1f} (i = {:.1f})'.format(lage_str, self.feh, self.rot, self.gdark_i*(180./np.pi))

        age_ind = self.age_index(lage)
        
        self.x = self.isocmds[age_ind][x_name]
        self.y = self.isocmds[age_ind][y_name]
        # For e.g. colors; assumes y is blue mag if CMD 
        if x_to_ymx:
            self.x = self.y - self.x
            # make x_name 'y_name - x_name'
            x_name = "{:s} - {:s}".format(y_name, x_name)

        self.x_name = x_name
        self.y_name = y_name

        if ax is not None:
            ax.set_ylabel(y_name)
            ax.set_xlabel(x_name)
            ax.set_xlim([x.max()+0.05, x.min()-0.05])
            ax.set_ylim([y.max()+0.2, y.min()-0.2])
            ax.set_title('MIST Isochrones: {:s} vs. {:s}'.format(y_name, x_name))

            # If true, this will plot a Geneva model at the given age for comparison, if 
            # the file for it exists.
            #if geneva_on:
                #plot_geneva_iso(ages[i], rots[i], pltmass, i, color_n, ax = ax)

        # stellar masses:
        self.init_masses = self.isocmds[age_ind]['initial_mass']

        # For masking out phases (e.g. PMS):
        if len(phasemask) > 0:
            self.phases = self.isocmds[age_ind]['phase']
            for i_pmask, pmask in enumerate(phasemask):
                unmasked_ind = np.where(phases != pmask)
                # Want a new array of only the valid phases on ea. pass to
                # keep indices matching between x, y, and p
                self.phases = phases[unmasked_ind]                    
                self.x = self.x[unmasked_ind]
                self.y = self.y[unmasked_ind]
                self.init_masses = self.init_masses[unmasked_ind]

        return #x, y, init_masses, phases

    def pltmass(self, ax, masses=None, alpha=1.0):

        """
            For plotting masses on an isochrone. Perhaps too complicated -- maybe delete in favor of something simpler.
        """
        x = self.x
        y = self.y
        init_masses = self.init_masses

        # pltmass may be given as a list of floats, or a float in order to plot a specific mass.
        if isinstance(masses, float):
            # Get index of nearest mass:
            diff_arr = abs(init_masses - masses)
            m_i = np.where(diff_arr == np.min(diff_arr))[0][0]
            ax.scatter(x[m_i], y[m_i], lw=0.1, alpha=0.5, color = base_line.get_color(), zorder=2)
            ax.text(x[m_i]*0.98, y[m_i], '{:.2f}'.format(init_masses[m_i]) + r' $M_{\odot}$', fontsize=4, color = base_line.get_color())

        # same as above, but for the case of a list of masses:
        elif isinstance(masses, list):
            for mass in pltmass:
                diff_arr = abs(init_masses - masses)
                m_i = np.where(diff_arr == np.min(diff_arr))[0][0]
                ax.scatter(x[m_i], y[m_i], lw=0.1, alpha=0.5, color = base_line.get_color(), zorder=2)
                ax.text(x[m_i]*0.98, y[m_i], '{:.2f}'.format(init_masses[m_i]) + r' $M_{\odot}$', fontsize=4, color = base_line.get_color())

        # If no masses are specified, plot some representative masses:
        else:
            base_line, = ax.plot(x, y, label=self.lbl, alpha = alpha)
            sc = ax.scatter(x, y, c = init_masses, s=6, lw=0.1, edgecolor=base_line.get_color(), cmap = plt.cm.gist_rainbow, alpha=0.5, zorder=-1)

            # Tracking individual masses for comparison of evolutionary states:
            tracked_masses = []
            massmatches = 0.0

            plotted_masses = []

            # The follwing for loop handles plotting several representative masses along the isochrone & labeling them
            int_masses = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
            for m_i, mass in enumerate(init_masses):

                # Only plot the full range of representative masses if on the first isochrone.
                near_int = (0.01 > np.amin(abs(int_masses - float('{:.2f}'.format(mass)))))
                # condition of brightest, bluest, or last point:
                brightest_bluest_last = ((x[m_i] == min(x)) or (y[m_i] == min(y)) or (m_i == len(init_masses) - 1))

                if  (near_int or brightest_bluest_last) & (i == 0):

                    # Labeling coordinates -- attempting to not overlap the isochrone & other labels.                    
                    txt_x = x[m_i] - 0.2
                    txt_y = y[m_i]
                    if m_i == len(init_masses) - 1:
                        txt_x += 0.25
                    elif y[m_i] == min(y):
                        txt_y -= 0.25
                        txt_x += 0.25

                    # text label for the mass value:
                    ax.text(txt_x, txt_y, '{:.2f}'.format(mass) + r' $M_{\odot}$', fontsize=4, color = base_line.get_color())
                    # Plots a scatter point at the location of the star w/ a particular mass:
                    ax.scatter(x[m_i], y[m_i], lw=0.1, alpha=0.5, color = base_line.get_color(), zorder=2)
                    # Store the masses used; plot these masses on subsequent isochrones:
                    tracked_masses.append(float('{:.2f}'.format(mass)))
            
                # To avoid crowding things, this segment plots mass labels on the side of the plot, rather than along the isochrone,
                # this occurs for all isochrones beyond the first one.
                # The masses picked out for plotting on the first isochrone are tracked and plotted if existing on other isochrones.
                nearint = (0.01 > np.amin(abs(int_masses - float('{:.2f}'.format(mass)))))
                notplt = float('{:.2f}'.format(mass)) not in plotted_masses
                tracked = (float('{:.2f}'.format(mass)) in np.unique(tracked_masses))

                nearint_notplt = (nearint) & (notplt)
                notplt_tracked = (notplt) & (tracked)

                if (nearint_notplt | notplt_tracked ) & (i > 0):

                    # massmatches increases y position to avoid overlap as more masses are found:
                    massmatches += 0.2
                    txt_x = min(x) - 0.3
                    txt_y = min(y) + massmatches + i*1.5

                    ax.text(txt_x-0.15, txt_y, '{:.2f}'.format(mass) + r' $M_{\odot}$', fontsize=4, color = base_line.get_color())
                    ax.scatter(x[m_i], y[m_i], lw=0.1, alpha=0.5, color = base_line.get_color(), zorder=2)

                    # Draw a line connecting the text to the data point:
                    ax.plot([txt_x, x[m_i]], [txt_y, y[m_i]], color=base_line.get_color(), lw=0.2, alpha=0.5)
                    # track plotted masses:
                    plotted_masses.append(float('{:.2f}'.format(mass)))

        return base_line, sc

    def isoplot(self, ax, masses=None, xlim=[], ylim=[], alpha=1.0, shade=0.0):
        # Plot a CMD of red vs. blue - red:
        if  shade > 1.0 or shade < 0.0:
            shade = 0.0
            print("Warning: shade must be a float between 0 and 1.")

        lc = plt.cm.Dark2(shade)

        base_line, = ax.plot(self.x, self.y, label=self.lbl, lw = 1, c = lc, alpha=alpha)
        if isinstance(masses, float):
            # Get index of nearest mass:
            diff_arr = abs(self.init_masses - masses)
            m_i = np.where(diff_arr == np.min(diff_arr))[0][0]
            ax.scatter(self.x[m_i], self.y[m_i], lw=0.1, alpha=0.5, color = base_line.get_color(), zorder=2)
            ax.text(self.x[m_i]*0.98, self.y[m_i], '{:.2f}'.format(self.init_masses[m_i]) + r' $M_{\odot}$', fontsize=8, color = 'k')
        
        #if self.xextent == [0,0]:
        #    self.xextent = [self.x.min()-0.05, self.x.max()+0.05]
        #if self.yextent == [0,0]:
        #    self.yextent = [self.y.max()+0.2, self.y.min()-0.2]

        if len(xlim) == 0:
            if self.xextent[1] < self.x.max()+0.05:
                self.xextent[1] = self.x.max()+0.05
            if self.xextent[0] > self.x.min()-0.05:
                self.xextent[0] = self.x.min()-0.05
        elif len(xlim) == 2:
            self.xextent = xlim

        if len(ylim) == 0: 
            if self.yextent[1] > self.y.min()-0.2:
                self.yextent[1] = self.y.min()-0.2 
            if self.yextent[0] < self.y.max()+0.2:
                self.yextent[0] = self.y.max()+0.2
        elif len(ylim) == 2:
            self.yextent = ylim

        ax.set_ylabel(self.y_name)
        ax.set_xlabel(self.x_name)
        ax.set_xlim(self.xextent)
        ax.set_ylim(self.yextent)
        ax.set_title('MIST Isochrones: {:s} vs. {:s}'.format(self.y_name, self.x_name))

        ax.legend(loc='best')

        return

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class EEP:
    
    """
    
    Reads in and plots MESA EEP files.

    
    """
    
    def __init__(self, filename, verbose=True):
        
        """
        
        Args:
            filename: the name of .track.eep file.
        
        Usage:
            >> eep = read_mist_models.EEP('00200M.track.eep')
            >> logTeff, center_h1, mdot = eep.eeps['log_Teff'], eep['center_h1'], eep['star_mdot']
            
        Attributes:
            version         Dictionary containing the MIST and MESA version numbers.
            abun            Dictionary containing Yinit, Zinit, [Fe/H], and [a/Fe] values.
            rot             Rotation in units of surface v/v_crit.
            minit           Initial mass in solar masses.
            hdr_list        List of column headers.
            eeps            Data.
            
        """
                        
        self.filename = filename

        if verbose:
            print('Reading in: ' + self.filename)
                        
        self.version, self.abun, self.rot, self.minit, self.hdr_list, self.eeps = self.read_eep_file()
        
    def read_eep_file(self):
        
        """

        Reads in the EEP file.
        
        Args:
            filename: the name of .track.eep file.
                
        """
        
        eeps = np.genfromtxt(self.filename, skip_header=11, names=True)
        
        with open(self.filename) as f:
            content = [line.split() for line in f]

        version = {'MIST': content[0][-1], 'MESA': content[1][-1]}
        abun = {content[3][i]:float(content[4][i]) for i in range(1,5)}
        rot = float(content[4][-1])
        minit = float(content[7][1])
        hdr_list = content[11][1:]
        
        return version, abun, rot, minit, hdr_list, eeps
        		
    def plot_HR(self, fignum=0, phases=[], phasecolor=[], phasemask = [], ax = None, cbar_valname=None, **kwargs):
        
        """

        Plots the HR diagram.

        Args:
            None.
            
        Keywords:
            accepts matplotlib keywords: color, linestyle, linewidth, etc.
            keyword: fignum, phase*, phasecolor, phasemask, ax, cbar_valname
            
            * Following the FSPS notation,
            * PMS:-1 ; MS:0 ; SGB+RGB:2 ; CHeB:3 ; EAGB:4 ; TPAGB:5 ; post-AGB:6 ; WR:9
    
        Usage:
            >> eep.plot_HR(fignum=3)
            >> eep.plot_HR(phase=[0, 2], phasecolor=['Gray', 'Blue']) #highlight the MS and RGB phases in gray and blue.

            >> eep.plot_HR(phasemask=[-1, 9], ax=ax1) # Do not plot the PMS or WR phase, and plot on existing Axes object 'ax1'. 
            >> eep.plot_HR(cbar_valname='surf_avg_omega_div_omega_crit') # Color lines by omega/omega_crit.
        
        """
        
        x = self.eeps['log_Teff']
        y = self.eeps['log_L']
        if isinstance(cbar_valname, str):
            cbar_vals = self.eeps[cbar_valname]
            if cbar_valname == 'surf_avg_omega_div_omega_crit':
                cbarmin = 0.0
                cbarmax = 1.0
            # placeholder...
            else:
                cbarmin = 0.0
                cbarmax = 1.0

            z = cbar_vals

        # For masking out phases (e.g. PMS):
        if len(phasemask) >= 0:
            p = self.eeps['phase']
            for i_pmask, pmask in enumerate(phasemask):
                unmasked_ind = np.where(p != pmask)
                # Want a new array of only the valid phases on ea. pass to
                # keep indices matching between x, y, and p
                p = p[unmasked_ind]                    
                x = x[unmasked_ind]
                y = y[unmasked_ind]
                if isinstance(cbar_valname, str):
                    z = z[unmasked_ind]
        
        if ax == None:
            fig = plt.figure(fignum)        
            ax = fig.add_subplot(111)

        ax.set_xlabel('log(Teff) [K]')
        ax.set_ylabel('log(L/Lsun)')


        # If supplied, create a color the lines by a third value:
        if isinstance(cbar_valname, str):
            ax, lc = color_lineby_z(x, y, z, cbarmin, cbarmax, ax=ax, alpha=kwargs['alpha'])
            ax.axis('auto')
            #ax.axis([max(x)+0.2, min(x)-0.2, min(y)-0.2, max(y)+0.2])

        # Or else just plot, don't color by 3rd value:
        else:
            ax.plot(x, y, **kwargs)
            #ax.axis([max(x)+0.2, min(x)-0.2, min(y)-0.2, max(y)+0.2])

        if len(phases) >= 0:

            phase_names = ['MS',None, 'SGB+RGB', 'CHeB', 'EAGB', 'TPAGB', 'post-AGB', 'WR', 'PMS']
            if len(phases) != len(phasecolor):
                print('The length of the phase and phasecolor array must be identical.')
                return
            for i_p, phase in enumerate(phases):
                # p may have been defined already when masking specified phases above, create it if not:
                if 'p' not in locals():
                    p = self.eeps['phase']

                p_ind = np.where(p == phase)
                if len(p_ind) > 0:
                    if phasecolor == '':
                        ax.plot(x[p_ind], y[p_ind], linewidth=4.0, alpha=0.5)
                        init_phase_age = self.eeps['star_age'][p_ind][0]
                        fin_phase_age = self.eeps['star_age'][p_ind][-1]
                        phase_deltaT = fin_phase_age - init_phase_age
                        #ax.text(x[p_ind][-1], y[p_ind][-1], r'{:s} $\delta$t = {:.2e}'.format(phase_names[phase], phase_deltaT), fontsize=3)
                    else:
                        ax.plot(x[p_ind], y[p_ind], color=phasecolor[i_p], linewidth=4.0, alpha=0.5)
                        init_phase_age = self.eeps['star_age'][p_ind][0]
                        fin_phase_age = self.eeps['star_age'][p_ind][-1]
                        phase_deltaT = fin_phase_age - init_phase_age
                        #ax.text(x[p_ind][-1], y[p_ind][-1], r'{:s} $\delta$t = {:.2e}'.format(phase_names[phase], phase_deltaT), fontsize=3)

        if isinstance(cbar_valname, str):
            return ax, lc, x, y
        else:
            return ax, x, y

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class EEPCMD:
    
    """
    
    Reads in and plots MESA EEP CMD files.

    
    """
    
    def __init__(self, filename, verbose=True):
        
        """
        
        Args:
            filename: the name of .track.eep.cmd file.
        
        Usage:
            >> eepcmd = read_mist_models.EEPCMD('00200M.track.eep.cmd')
            >> B, V, mdot = eepcmd.eepcmds['Bessell_B'], eep['Bessell_V'], eep['star_mdot']
            
        Attributes:
            version         Dictionary containing the MIST and MESA version numbers.
            photo_sys       Photometric system.
            abun            Dictionary containing Yinit, Zinit, [Fe/H], and [a/Fe] values.
            rot             Rotation in units of surface v/v_crit.
            minit           Initial mass in solar masses.
            hdr_list        List of column headers.
            Av_extinction   Av for CCM89 extinction.
            eepcmds         Data.
            
        """
                        
        self.filename = filename
        if verbose:
            print('Reading in: ' + self.filename)
                        
        self.version, self.photo_sys, self.abun, self.rot, self.minit, self.Av_extinction, self.hdr_list, self.eepcmds = self.read_eepcmd_file()
        
    def read_eepcmd_file(self):
        
        """

        Reads in the EEP CMD file.
        
        Args:
            filename: the name of .eep.cmd file.
                
        """
        
        eepcmds = np.genfromtxt(self.filename, skip_header=14, names=True)
        
        with open(self.filename) as f:
            content = [line.split() for line in f]

        version = {'MIST': content[0][-1], 'MESA': content[1][-1]}
        photo_sys = ' '.join(content[2][4:])
        abun = {content[4][i]:float(content[5][i]) for i in range(1,5)}
        rot = float(content[5][-1])
        minit = float(content[8][1])
        Av_extinction = float(content[11][-1])
        hdr_list = content[14][1:]
        
        return version, photo_sys, abun, rot, minit, Av_extinction, hdr_list, eepcmds
        		
    def plot_CMD(self, filters, fignum=0, phases=[], phasecolor=[], **kwargs):
        
        """

        Plots the CMD diagram.

        Args:
            filters: a list of three filters, ['filter1', 'filter2', 'filter3']. x-axis: 'filter1'-'filter2', y-axis: 'filter3'
            
        Keywords:
            accepts matplotlib keywords: color, linestyle, linewidth, etc.
            keyword: fignum, phase*, phasecolor
            
            * Following the FSPS notation,
            * PMS:-1 ; MS:0 ; SGB+RGB:2 ; CHeB:3 ; EAGB:4 ; TPAGB:5 ; post-AGB:6 ; WR:9
    
        Usage:
            >> eepcmd.plot_CMD(['Bessell_B', 'Bessell_V', 'Bessell_V'], fignum=3)
        """
        
        try:
            x1 = self.eepcmds[filters[0]]
        except:
            print(filters[0]) + ' does not appear in this file.'
            return
        try:
            x2 = self.eepcmds[filters[1]]
        except:
            print(filters[1]) + ' does not appear in this file.'
            return
        try:
            y = self.eepcmds[filters[2]]
        except:
            print(filters[2]) + ' does not appear in this file.'
            return
        
        fig = plt.figure(fignum)
        plt.xlabel(' '.join(filters[0].split('_')) + '-' + ' '.join(filters[1].split('_')), fontsize=22)
        plt.ylabel(' '.join(filters[2].split('_')), fontsize=22)
        
        ax = fig.add_subplot(111)
        ax.plot(x1-x2, y, **kwargs)
        ax.axis([min(x1-x2)-0.2, max(x1-x2)+0.2, max(y)+0.2, min(y)-0.2])

        if len(phases) >= 0:
            if len(phases) != len(phasecolor):
                print('The length of the phase and phasecolor array must be identical.')
                return
            for i_p, phase in enumerate(phases):
                p = self.eepcmds['phase']
                p_ind = np.where(p == phase)
                if len(p_ind) > 0:
                    if phasecolor == '':
                        ax.plot(x1[p_ind]-x2[p_ind], y[p_ind], linewidth=4.0, alpha=0.5)
                    else:
                        ax.plot(x1[p_ind]-x2[p_ind], y[p_ind], color=phasecolor[i_p], linewidth=4.0, alpha=0.5)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def color_lineby_z(x, y, z, zmin, zmax, ax, alpha=1.0):
        
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, cmap=plt.get_cmap('jet'))
    lc.set_array(z)
    lc.set_clim(vmin=0.0,vmax=1.0)
    lc.set_alpha(alpha)

    ax.add_collection(lc)

    return ax, lc
