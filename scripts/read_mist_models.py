from __future__ import print_function
import numpy as np
from matplotlib.collections import LineCollection
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import glob
import seaborn as sns
from collections import OrderedDict
from pltstyle.plotstyle import set_style
set_style()
# no seaborn gridlines on plot, and white bg:
#sns.set_context('paper')
#sns.set(font='serif')
#sns.set_style("ticks", {'font.family' : 'serif', 'font.serif': ['Times', 'Palatino', 'serif']})
#plt.axes(frameon=False)

from MIST_codes.scripts.fileio import *
from MIST_codes.scripts.gdiso import *
from MIST_codes.scripts import read_geneva as rg
from MIST_codes.scripts import isomist

#params = {'axes.labelsize': 20,'axes.titlesize':20, 'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14}
#mpl.rcParams.update(params)

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText
from matplotlib.offsetbox import AnchoredText
from matplotlib.patches import Rectangle

# plot font sizing:
mpl.rcParams.update(set_style())
#{'axes.labelsize':20, 'axes.titlesize': 20, 'font.size':14,
          #'xtick.labelsize': 14, 'ytick.labelsize': 14, 'font.family':"serif", 'mathtext.fontset':'stix'}
#mpl.rcParams.update(params)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ISO:
    
    """
    
    Reads in MIST isochrone files.

    
    """
    
    def __init__(self, feh = 0.00, vvcrit = 0.0, filename = None, gravdark_i = 0.0, exttag = None, verbose=False, read=True, version='1.0'):
    
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
        
        if filename != None:
            self.filename = filename
        else:
            filename = get_fn(feh, vvcrit, mode='iso', gravdark_i = gravdark_i, exttag = exttag, version = version)
            #print(filename)
            self.filename = filename

        if verbose and read:
            print('Reading in: ' + self.filename)
        if read:
            self.version, self.abun, self.rot, self.ages, self.num_ages, self.hdr_list, self.isos = self.read_iso_file()
        # Check if gravity darkening is desired:
        if gravdark_i > 0.0:
            #print("Checking for gravity darkened .iso file...")
            # Will check for GDed .iso file; creates one if nec.
            gdisoname = filename.split('.iso')[0] + '_gdark{:.1f}.iso'.format(gravdark_i)
            # Check if the .iso file already exists:
            if os.path.isfile(gdisoname):
                #print("Gravity darkened .iso file {:s} exists.".format(gdisoname))
                self.filename = gdisoname#outfname
                if read:
                    self.version, self.abun, self.rot, self.ages, self.num_ages, self.hdr_list, self.isos = self.read_iso_file()
            # if it doesn't, create it (have to read in isos to calculate gd'ed values):
            else:
                #print("Creating gravity darkened .iso file {:s}.".format(gdisoname))
                self.version, self.abun, self.rot, self.ages, self.num_ages, self.hdr_list, self.isos = self.read_iso_file()
                #print("Loaded non-gravity darkened .iso file {:s}. Re-writing with gravity darkened values.".format(filename))
                filename = rw_gdiso(self, gdisoname, gravdark_i)
                #print("Created gravity darkened .iso file {:s}.".format(gdisoname))             
                self.filename = filename
                
            # replace filename with the new GD file:
            #self.filename = filename
            # read in the new gravity darkened data:
            #if read:
            #    self.version, self.abun, self.rot, self.ages, self.num_ages, self.hdr_list, self.isos = self.read_iso_file()

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
        print(len(iso_set))
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
    
    def __init__(self, feh=0.00, vvcrit=0.0, ebv=0.0, gravdark_i=0.0, exttag=None, filename=None, verbose=False, version = '1.0'):
    
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
            ages            List of unique log10 ages.
            num_ages        Number of ages.
            hdr_list        List of column headers.
            isocmds         Data (data columns corresp. to column headers).
        
        """
        Av = 3.1*ebv
        if filename != None:
            if '.cmd' not in filename:
                iso_filename = ISO(feh = feh, vvcrit = vvcrit, gravdark_i = gravdark_i, exttag = exttag, read=False, filename=filename, version = version).filename
                self.filename = isomist.createcmd(iso_filename, Av = Av, gravdark_i = gravdark_i)
            else:
                self.filename = filename
        else:
            iso_filename = ISO(feh = feh, vvcrit = vvcrit, gravdark_i = gravdark_i, exttag = exttag, read=False, filename=filename, version = version).filename
            self.filename = isomist.createcmd(iso_filename, Av = Av, gravdark_i = gravdark_i)


        self.exttag = exttag
        self.feh = feh
        self.gdark_i = gravdark_i
        self.xextent = [0,0]
        self.yextent = [0,0]
 
        if verbose:
            print('Reading in: ' + self.filename)
        
        self.version, self.photo_sys, self.abun, self.Av_extinction, self.rot, self.ages, self.num_ages, self.hdr_list, self.isocmds = self.read_isocmd_file()
        #print("Created an ISOCMD instance with photoemetric system {:s}, [Fe/H] = {:f}, Av = {:f}, V/Vc = {:.1f}, inclination angle = {:.1f}".format(self.photo_sys, self.feh, self.Av_extinction, self.rot, self.gdark_i))
    
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

        Returns the index for the user-specified log10 age.

        usage:
        
            isocmd.age_index(8.00)
        
        Args:
            age: the log10 age of the isochrone.
        
        """
        
        diff_arr = abs(np.array(self.ages) - age)
        age_index = np.where(diff_arr == min(diff_arr))[0][0]
        
        if ((age > max(self.ages)) | (age < min(self.ages))):
            print('The requested age is outside the range. Try between ' + str(min(self.ages)) + ' and ' + str(max(self.ages)))
            
        return age_index

    def get_masked(self, hdr_name, phasemask, age_ind=slice(None,None), dmod=0.0, mrange=[]):

        # For masking out phases (e.g. PMS); grabs appropriate age block if told to:
        x = self.isocmds[age_ind][hdr_name]
        p = self.isocmds[age_ind]['phase']
        miraw = self.isocmds[age_ind]['initial_mass']
        mi = np.array(list(map(float, miraw)))

        # cut data to a given mass range:
        if len(mrange) == 2:
            mrange_i = (min(mrange) <= mi) & (mi <= max(mrange))
            x = x[mrange_i]
            p = p[mrange_i] 

        # if phasemask = [] (default), this loop is skipped and there's no mask.
        for pmask in phasemask:
            unmasked_ind = np.where(p != pmask)
            # Want a new array of only the valid phases on ea. pass to
            # keep indices matching between x, y, and p
            p = p[unmasked_ind]                    
            x = x[unmasked_ind]

        masked_data = x + dmod

        return masked_data

    def get_data(self, hdr_names, phasemask, lage=None, age_ind=None, dmod=0.0, mrange=[]):

        # returns a dictionary of data with keys as converted header names and items as the corresp. data columns.

        if age_ind == None:
            if lage != None:
                age_ind = self.age_index(lage)
            else:
                print('To gather isochrone data, please specify the isochrones log10 age (paramater lage not set).')
                print('Alternatively, please supply the data array index corresponding to the desired age block.')
                return
        

        # get data values for given x, y names (considering phase masks).
        data_list = []
        data_dict = {}
        for name in hdr_names:

            # check for math operations (all of these get the data for the various cases); only works for ops btwn two columns atm.
            # e.g. you give hdr_names = ['Bessel_B - Bessel_V', 'Bessel_V'], this will set a data element as B-V with key 'B-V' and 
            # another w/ V and key 'V'.
            # will mask evolutionary phases accoring to 'phasemask' argument.
            if '-' in name:
                splitnames = name.split('-')
                data = self.get_masked(splitnames[0], phasemask=phasemask, age_ind=age_ind, dmod=dmod, mrange=mrange) - \
                       self.get_masked(splitnames[1], phasemask=phasemask, age_ind=age_ind, dmod=dmod, mrange=mrange)
            elif '+' in name:
                splitnames = name.split('+')
                data = self.get_masked(splitnames[0], phasemask=phasemask, age_ind=age_ind, dmod=dmod, mrange=mrange) + \
                       self.get_masked(splitnames[1], phasemask=phasemask, age_ind=age_ind, dmod=dmod, mrange=mrange)
            elif '*' in name:
                splitnames = name.split('*')
                data = self.get_masked(splitnames[0], phasemask=phasemask, age_ind=age_ind, dmod=dmod, mrange=mrange) * \
                       self.get_masked(splitnames[1], phasemask=phasemask, age_ind=age_ind, dmod=dmod, mrange=mrange)
            elif '/' in name:
                splitnames = name.split('/')
                data = self.get_masked(splitnames[0], phasemask=phasemask, age_ind=age_ind, dmod=dmod, mrange=mrange) / \
                       self.get_masked(splitnames[1], phasemask=phasemask, age_ind=age_ind, dmod=dmod, mrange=mrange)
            else:
                data = self.get_masked(name, phasemask=phasemask, age_ind=age_ind, dmod=dmod, mrange=mrange)

            # will store the data column in a dictionary with the key as the column's corresponding name (converted).
            data_list.append((convert_name(name), data))

        # OrderedDict remembers the order in which items are placed in the dictionary. Uses a list of tuples corresp. to the dictionary's items.
        data_dict = OrderedDict(data_list)

        return data_dict

    def set_isodata(self, lage, x_name, y_name, dmod=0.0, ax=None, phasemask=[], 
                    x_to_ymx=True, geneva_on=False, labels=['exttag', 'age', 'feh', 'vvc', 'inc'], 
                    texts=None, lgdlbl=None, verbose=False):
        """
            Given a log10 age, get the data desired to plot the CMD or HRD
           
            Control labeling of isochrones:

                Use the "labels" argument; this should be a list, valid elements are:
                    + exttag, age, feh, vvc, inc
                Equate this argument to a list containing either of the above elements as a string, e.g.
                    labels = ['feh']
                Neglect any one or more of the elements to exclude it from the label.
            
        """

        exttag = self.exttag

        age_ind = self.age_index(lage)

        if verbose:
            print("Given log10 age: {:.2f}; Closest log10 age: {:.2f}".format(lage, self.isocmds[age_ind]['log10_isochrone_age_yr'][0]))
            print("Using log10 age = {:.2f}.".format(self.isocmds[age_ind]['log10_isochrone_age_yr'][0]))
        lage = self.isocmds[age_ind]['log10_isochrone_age_yr'][0]

        # For labeling age in potential plots, :
        if lage >= 9.0:
            lage_str = "{:.1f} Gyr".format((10**lage)/(10.**9))
        if 6.0 <= lage < 9.0:
            lage_str = "{:.1f} Myr".format((10**lage)/(10.**6))
  
        # extra tag conversions (to less eye soring versions):
        if 'mlta' in exttag:
            exttag = r"$\alpha_{\rm{MLT}}=$"+"{:s}".format(exttag.split('mlta')[-1]) 
   
        # Creates a plot (legend) label for the isochrone.
        if lgdlbl == None:
            lbl = 'MIST({:s}): age = {:s}, [Fe/H] = {:.2f}, '.format(exttag, lage_str, self.feh) + \
                          r'$\frac{\Omega}{\Omega_c}$' + '  = {:.1f} (i = {:.1f})'.format(self.rot, self.gdark_i)
            if 'exttag' not in labels and exttag != None:
                lbl = lbl.replace('({:s})'.format(exttag), '')
            if 'age' not in labels:
                lbl = lbl.replace('age = {:s}, '.format(lage_str), '')
            if 'feh' not in labels:
                lbl = lbl.replace('[Fe/H] = {:.2f}, '.format(self.feh), '')
            if 'vvc' not in labels:
                lbl = lbl.replace(r'$\frac{\Omega}{\Omega_c}$' + ' = {:.1f} '.format(self.rot), '')
            if 'inc' not in labels:
                lbl = lbl.replace('(i = {:.1f})'.format(self.gdark_i), '')
            # check for cleanup of leftover comma and colon:
            lbl = lbl.replace(', (', ' (')
            lbl = lbl.replace(': (', ' (')
            if len(labels) == 1 and labels[0] == 'exttag':
                lbl = lbl.replace(':', '')
            self.lbl = lbl
        else:
            self.lbl = lgdlbl
        
        # for creating text labels:
        if texts is not None:
            if isinstance(texts, list):
                #textlbl = ''
                strs = []
                for ele in texts:
                    if ele == 'age':
                        strs.append('age = {:s}'.format(lage_str))
                    elif ele =='feh' in texts:
                        strs.append('[Fe/H] = {:.2f}'.format(self.feh))
                    elif ele == 'vvc' in texts:
                        strs.append(r'$\frac{\Omega}{\Omega_c}$' + ' = {:.1f}'.format(self.rot))
                    elif ele == 'inc' in texts:
                        strs.append('i = {:.1f}'.format(self.gdark_i)+r'$^{\circ}$')
                    else:
                        strs.append('{:s}'.format(ele))
                self.txtlbl = strs

#                else:
#                    for ele in texts:
#                        if ele == 'age':
#                            textlbl += 'age = {:s}\n'.format(lage_str)
#                        elif ele =='feh' in texts:
#                            textlbl += '[Fe/H] = {:.2f}\n'.format(self.feh)
#                        elif ele == 'vvc' in texts:
#                            textlbl += r'$\frac{\Omega}{\Omega_c}$' + '  = {:.1f}\n'.format(self.rot)
#                        elif ele == 'inc' in texts:
#                            textlbl += 'i = {:.1f}\n'.format(self.gdark_i)
#                        else:
#                            textlbl += '{:s}\n'.format(ele)
#
#                self.txtlbl = textlbl
            
            # if a ready made string is provided:
            elif isinstance(texts, str):
                self.txtlbl = [texts]
        else:
            self.txtlbl = None 

        # datadict holds all .iso file data for the specified isochrone. Organize it into desired x, y here:
        datdict = self.get_data([x_name, y_name], phasemask, age_ind=age_ind, dmod=dmod)
        self.x_name, self.y_name = datdict.keys()
        self.x, self.y = list(datdict.values())

        # get initial masses
        self.init_masses = list((self.get_data(['initial_mass'], phasemask, age_ind=age_ind)).values())[0]
        self.init_masses = np.array(list(map(float, self.init_masses)))

        # sets labels if given an axis...not necessary?:
        if ax is not None:
            ax.set_ylabel(u"${:s}$".format(self.y_name), rotation=90)
            ax.set_xlabel(u"${:s}$".format(self.x_name))
            ax.set_xlim([x.max()+0.05, x.min()-0.05])
            ax.set_ylim([y.max()+0.2, y.min()-0.2])

        # stellar masses (not needed?):
        # self.init_masses = self.isocmds[age_ind]['initial_mass']

        if verbose:
            print("Set data for {:s}".format(self.lbl))

        return #x, y, init_masses, phases

    def pltmass(self, ax, masses=None, alpha=1.0):

        """
            For plotting masses on an isochrone. Perhaps too complicated -- maybe delete in favor of something simpler.
            
            *** DEPRECATED? ***

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
            if not init_masses[m_i].is_integer():
                ax.text(x[m_i]*0.98, y[m_i], '{:.2f}'.format(init_masses[m_i]) + r' $M_{\odot}$', fontsize=4, color = base_line.get_color())
            else:
                ax.text(x[m_i]*0.98, y[m_i], '{:d}'.format(int(init_masses[m_i])) + r' $M_{\odot}$', fontsize=4, color = base_line.get_color())

        # same as above, but for the case of a list of masses:
        elif isinstance(masses, list):
            for mass in pltmass:
                diff_arr = abs(init_masses - masses)
                m_i = np.where(diff_arr == np.min(diff_arr))[0][0]
                ax.scatter(x[m_i], y[m_i], lw=0.1, alpha=0.5, color = base_line.get_color(), zorder=2)
                ax.text(x[m_i]*0.98, y[m_i], '{:.2f}'.format(init_masses[m_i]) + r' $M_{\odot}$', fontsize=4, color = base_line.get_color())

        return base_line, sc

    def isoplot(self, ax, masses=None, xlims=None, ylims=None, shade=None, legend=False, label=False, setlim=False, 
                title=False, obs=None, justify_lbl=None, lcon=False, **kwargs):

        """
            Plots CMD of a MIST isochrone object.

            Args:
            -----
                + ax: Matplotlib.pyplot axes object -- pass an axes and the plot goes on it.
                + masses: Optionally specify masses to mark on the isochrone plot.
                + xlim: Give a list of two numbers to specify x axis limits.
                + ylim: Same as xlim but for y axis.
                + shade: A number between 0 and 1 to choose a line color -- goes through some cmap used by 'lc' below.
                + legloc: legend location -- uses matplotlib legend() method keywords.
                + legsize: Size of the legend.
                + **kwargs: any keyword arguments for matplotlib.pyplot.plot() to be used when plotting the isochrone.
        """

        # Plot a CMD of red vs. blue - red:
        # Alter shade of line:
        if shade != None:
            if  shade > 1.0 or shade < 0.0:
                shade = 0.0
                print("Warning: shade should be between 0 and 1; defaulting to 0.0.")

            lc = plt.cm.Dark2(shade)
            kwargs['c'] = lc

        # Add (predefined) text label (deprecate?):
        if isinstance(label, bool):
            if label == True:                
                # handles text labels via AnchoredText (self.txtlbl set in set_isodata()):
                if self.txtlbl != None:
                    #at = AnchoredText(self.txtlbl, **label)
                    proxarts.append(Rectangle((0,0), 1, 1, fc='w', fill=False, edgecolor='none', alpha=0.0, linewidth=0, label=self.txtlbl))
                    lg = ax.legend(handles=proxarts, **label)
                    # right/left align text
                    if justify_lbl == 'right':
                        vp = lg._legend_box._children[-1]._children[0]
                        vp.align ='right'
                    elif justify_lbl == 'left':
                        vp = lg._legend_box._children[-1]._children[0]
                        for c in vp._children:
                            c._children.reverse()
                        vp.align ='left'

                    ax.add_artist(lg)

        elif isinstance(label, dict):
            # label may be passed as a dictionary for AnchoredText kwargs; 
            # e.g. label={"frameon":False,"prop":dict(size=8),"loc":2}
            if self.txtlbl != None:
                proxarts = []
                if len(self.txtlbl) == 1:
                    lbls = (self.txtlbl[0]).split('\n')
                elif len(self.txtlbl) > 1:
                    lbls = self.txtlbl
                for lbl in lbls:
                    # invisible artists for text labels
                    proxarts.append(Rectangle((0,0), 1, 1, edgecolor='none',alpha=0.0, linewidth=0, label=lbl))
                if lcon:
                    proxarts.append(Rectangle((0,0), 1, 1, facecolor=kwargs['c'], edgecolor='none',alpha=0.6, linewidth=0, label=''))             

                lg = ax.legend(handles=proxarts, **label)
                # right/left align text
                if justify_lbl == 'right':
                    print('aligning right')
                    vp = lg._legend_box._children[-1]._children[0]
                    if lcon:
                        vp._children[-1]._children.reverse()
                    vp.align ='right'
                elif justify_lbl == 'left':
                    vp = lg._legend_box._children[-1]._children[0]
                    for c in vp._children:
                        c._children.reverse()
                    if lcon:
                        vp._children[-1]._children.reverse()
                    vp.align ='left'

                print(ax)
                # add legend w/ text labels to ax:
                ax.add_artist(lg)

        # line width = 1
        kwargs['lw'] = 1
        # legend label
        kwargs['label'] = self.lbl

        # Plot the CMD:
        
        #with data if a data file is specified:
        if obs is not None:
            mags = np.genfromtxt(os.path.join('/n/conroyfs1/sgossage/match_photometry/TIGS/phot_files', '{:s}.phot'.format(obs)))
            v, i = mags.T
            ax.scatter(v-i, i)

        base_line, = ax.plot(self.x, self.y, **kwargs)

        if setlim:
            ax.set_ylim([min(self.y), max(self.y)])
            ax.set_xlim([min(self.x), max(self.x)])

        # check if axis limits need to be expanded to fit new data; expand if so.
        curr_xlims, curr_ylims = expand_lims(ax, x=self.x, y=self.y, xlims=xlims, ylims=ylims, reverse_y=True)

        # Labeling masses if they are provided in a list.
        if isinstance(masses, float):
            # Get index of nearest mass:
            diff_arr = abs(self.init_masses - masses)
            m_i = np.where(diff_arr == np.min(diff_arr))[0][0]
            ax.scatter(self.x[m_i], self.y[m_i], lw=0.1, alpha=0.5, color = base_line.get_color(), zorder=2)
            print('!')
            if not round(self.init_masses[m_i], 2).is_integer():
                ax.text(self.x[m_i]*0.98, self.y[m_i], '{:.2f}'.format(self.init_masses[m_i]) + r' $M_{\odot}$', fontsize=8, color = 'k')
            else:
                ax.text(x[m_i]*0.98, y[m_i], '{:d}'.format(int(self.init_masses[m_i])) + r' $M_{\odot}$', fontsize=4, color = base_line.get_color())

            
        curr_xlims, curr_ylims = expand_lims(ax, ypad=0.05, xpad=0.05)

        self.xextent = np.array(curr_xlims)
        self.yextent = np.array(curr_ylims)

        ax.set_ylabel(u"${:s}$".format(self.y_name))
        ax.set_xlabel(u"${:s}$".format(self.x_name))
        if title:
            ax.set_title(u'MIST Isochrones: ${:s}$ vs. ${:s}$'.format(self.y_name, self.x_name))

        if isinstance(legend, bool):
           if legend:
               legend = ax.legend(loc='best', prop={'size':8}, frameon=False)
        elif isinstance(legend, dict):
            # legend may be passed as a dictionary for legend() kwargs; 
            # e.g. label={"frameon":False,"prop":dict(size=8),"loc":'upper right'}
            legend = ax.legend(**legend)

        return self.x, self.y, base_line

    def colormi(self, limits, ax=None, cmap=plt.cm.cool):

        mrange_i = (min(limits) <= self.init_masses) & (self.init_masses <= max(limits))

        mi = self.init_masses[mrange_i]
        print(mi)
        x = self.x[mrange_i]
        y = self.y[mrange_i]
        if ax is not None:
            ax, lc = color_lineby_z(x, y, mi, min(mi), max(mi), ax, cmap=cmap)

        return x, y
        

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class EEP:
    
    """
    
    Reads in and plots MESA EEP files.

    
    """
    
    def __init__(self, feh=0.00, vvcrit=0.0, mass=1.0, filename=None, gravdark_i = 0.0, exttag = None, verbose=True):
        
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
        
        # if a track.eep file name is given, use it.
        if filename != None:
            self.filename = filename

        # or else find a track.eep file given feh, v/vcrit, mass, and an extra tag (e.g. 'TP'):
        else:
            print('Getting file name...')
            self.filename = get_fn(feh, vvcrit, mode='eep', mass=mass, exttag=exttag)
            if isinstance(self.filename, list):
                # find the closest mass:
                masslist = self.filename
                masslist = list(map(float, masslist))
                diffarr = abs(np.array(masslist) - mass)
                close_idx = np.where(diffarr == diffarr.min())[0][0]
                mass = masslist[close_idx]
                print(mass)
                print('Using {:.2f} Msol instead.'.format(mass))
                self.filename = get_fn(feh, vvcrit, mode='eep', mass=mass, exttag=exttag)

        if verbose:
            print('Reading in: ' + self.filename)

        # dictionary for phases & their names:
        self.phase_names = {-1:'PMS', 0:'MS', 2:'SGB+RGB', 3:'CHeB', 4:'EAGB', 5:'TPAGB', 6:'post-AGB', 9:'WR'}
        
        # reads the MIST version number, abundances, v/vcrit, initial mass, header list, data (eeps):
        self.version, self.abun, self.rot, self.minit, self.hdr_list, self.eeps = self.read_eep_file()

        if gravdark_i > 0.0:
            # rewrite the luminosity and effective temp according to grav. darkening:
            self.filename = rw_gdeep(self, self.filename, gravdark_i)
            self.version, self.abun, self.rot, self.minit, self.hdr_list, self.eeps = self.read_eep_file()

        self.gravdark_i = gravdark_i

        self.lbl = r'MIST: M = {:.2f}'.format(self.minit)+r' $M_{\odot}$, '+r'$\Omega/\Omega_c$'+' = {:.1f}, i = {:.2f} deg'.format(self.rot, self.gravdark_i)

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

    def get_masked(self, hdr_names, phasemask):

        """
            Returns a data column that has had specified phases removed from it (i.e. removes rows corresp. to those undesired
        phases of evolution).
        """

        masked_data = []
        for name in hdr_names:
            # For masking out phases (e.g. PMS):
            x = self.eeps[name]
            p = self.eeps['phase']
            for pmask in phasemask:
                unmasked_ind = np.where(p != pmask)
                # Want a new array of only the valid phases on ea. pass to
                # keep indices matching between x, y, and p
                p = p[unmasked_ind]                    
                x = x[unmasked_ind]

            masked_data.append(x)

        return masked_data

    def get_phaselen(self, phasemask):
        """
            Returns the length of time of a phase of evolution in years.
        """
        ages = self.get_masked(['star_age'], phasemask=phasemask)[0]
        return ages[-1] - ages[0]
        		
    def plot_HR(self, fignum=0, phases=[], phasecolor=[], phasemask = [], ax = None, shade = 0.0, cbar_valname = None, masslbl=True, vvclbl=False,
                title=None, xlims = None, ylims = None, agemarks = None, legend = None, label=False, showplt = False, savename = None, xlabel=True,
                ylabel = True, geneva_on=False, geneva_lbl=False, MIST_on = True, setlim=False, **kwargs):
        
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

        if phasemask == 'allbutms':
            phasemask = [-1,2,3,4,5,6,9]
        self.phasemask = phasemask

        if shade != None:
            if  shade > 1.0 or shade < 0.0:
                shade = 0.0
                print("Warning: shade must be between 0 and 1; defaulting to 0.0.")

            lc = plt.cm.Dark2(shade)
            kwargs['c'] = lc

        if not MIST_on:
            kwargs['alpha'] = 0.0

        hdr_names = ['log_Teff', 'log_L']
        if cbar_valname is not None:
            hdr_names.append(cbar_valname)
        #x = self.eeps['log_Teff']
        #y = self.eeps['log_L']
        if isinstance(cbar_valname, str):
            cbar_vals = self.eeps[cbar_valname]
            if cbar_valname == 'surf_avg_omega_div_omega_crit':
                cbarmin = 0.0
                cbarmax = 1.0
            # placeholder...should adapt to cbar parameter value range?
            else:
                cbarmin = 0.0
                cbarmax = 1.0

        #    z = cbar_vals

        # For masking out phases (e.g. PMS):
        if len(phasemask) > 0:
            if isinstance(cbar_valname, str):
                x, y, z = self.get_masked(hdr_names, phasemask=phasemask) 
            else:
                x, y = self.get_masked(hdr_names, phasemask=phasemask)
            x = x[::-1]
            y = y[::-1]
        else:
            x = self.eeps['log_Teff'][::-1]
            y = self.eeps['log_L'][::-1]
        
        if ax == None:
            fig = plt.figure(fignum)        
            ax = fig.add_subplot(111)

        if xlabel:
            ax.set_xlabel(r'$\rm{log}(T_{eff})$')
        if ylabel:
            ax.set_ylabel(r'$\rm{log}(L/L_{\odot})$')
    
        if label == True:
            kwargs['label'] = self.lbl
        elif isinstance(label, str):
            print(label)
            kwargs['label'] = label

        # If supplied, create a color the lines by a third value:
        if isinstance(cbar_valname, str):
            ax, lc = color_lineby_z(x, y, z, cbarmin, cbarmax, ax=ax, **kwargs)#alpha=kwargs['alpha'])
            ax.axis('auto')
            #ax.axis([max(x)+0.2, min(x)-0.2, min(y)-0.2, max(y)+0.2])

        # Or else just plot, don't color by 3rd value:
        else:
            ax.plot(x, y, **kwargs)
            #ax.axis([max(x)+0.2, min(x)-0.2, min(y)-0.2, max(y)+0.2])

        if len(phases) >= 0:

            phase_names = self.phase_names
            if len(phases) != len(phasecolor):
                print('The length of the phase and phasecolor array must be identical.')
                return
            for i_p, phase in enumerate(phases):
                # p may have been defined already when masking specified phases above, create it if not:
                if 'p' not in locals():
                    p = self.eeps['phase']

                p_ind = np.where(p == phase)[0][::-1]
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

        self.x = x
        self.y = y

        # If desired, plot an analogous Geneva model for comparison:

        # maintain current axis limits prior to Geneva plotting.
        if setlim:
            ax.set_ylim([min(self.y), max(self.y)])
            ax.set_xlim([min(self.x), max(self.x)])

        curr_ylims = ax.get_ylim()
        curr_xlims = ax.get_xlim()

        if geneva_on:
            try:
                geneva_star = rg.star(vvc = self.rot, minit = self.minit)
                if not MIST_on:
                    try:
                        kwargs['alpha'] = alpha
                    except NameError:
                        kwargs['alpha'] = 1.0
                        
                    # use provided kwargs for geneva plot if mist is not turned on.
                    genx, geny, lc = geneva_star.plot_HR(ax = ax, label = geneva_lbl, shade = shade, ls = '--', **kwargs)
                else:
                    genx, geny, lc = geneva_star.plot_HR(ax = ax, label = geneva_lbl, shade = shade, ls = '--')

            # in case the Geneva star is not recoverable:
            except IOError:
                print('Failed to create Geneva star: v/vcrit = {:.1f}, Mass = {:.2f}'.format(self.rot, self.minit))

        ax.set_xlim(curr_xlims)
        ax.set_ylim(curr_ylims)
        curr_xlims, curr_ylims = expand_lims(ax, x=self.x, y=self.y, xlims=xlims, ylims=ylims, reverse_x=True)

        # mark the tracks if desired:
        if masslbl or vvclbl:
            curr_xlims, curr_ylims = self.mark_track(ax, xlims = xlims, ylims = ylims, masslbl=masslbl, vvclbl=vvclbl)
        
        # apply a 1% pad to the axes:
        curr_xlims, curr_ylims = expand_lims(ax, ypad=0.01, xpad=0.001)
     
        # marking age points on track:
        if agemarks != None:
            ages = self.eeps['star_age']
            oms = self.eeps['surf_avg_omega_div_omega_crit']
            # For masking out phases (e.g. PMS):
            if len(phasemask) >= 0:
                p = self.eeps['phase']
                for i_pmask, pmask in enumerate(phasemask):
                    unmasked_ind = np.where(p != pmask)
                    # Want a new array of only the valid phases on ea. pass to
                    # keep indices matching between x, y, and p
                    p = p[unmasked_ind]                    
                    #x = x[unmasked_ind]
                    #y = y[unmasked_ind]
                    ages = ages[unmasked_ind]
                    oms = oms[unmasked_ind]

                for agemark in agemarks:
                    agemark = 10**agemark
                    dage_arr = abs(ages - agemark)
                    agediff_i = np.where(dage_arr == np.amin(dage_arr))[0][0]
                    agemx = x[agediff_i]
                    agemy = y[agediff_i]
                    ax.scatter(agemx, agemy, marker = '+', s=20, c = 'k', lw = 0.8)
                    print("{:.1f} at {:.3f}: {:.2f}".format(self.rot, ages[agediff_i], oms[agediff_i]))
                        
            
        if isinstance(cbar_valname, str):
            cbar = plt.colorbar(lc, ax=ax)
            cbar.set_label(cbar_valname)
        #--------------------------------------------------------------------------------------

        if title != None:
            ax.set_title(title)

        if legend != None:
            if isinstance(legend, bool):
                ax.legend(loc = legloc, prop={'size':8}, frameon=False)
            elif isinstance(legend, dict):
                ax.legend(**legend)

        if savename != None:
            plt.savefig(savename, dpi=600)
        if showplt:
            print(ax.get_ylim())
            print(ax.get_xlim())
            plt.show()

        if isinstance(cbar_valname, str):
            return lc, x, y
        else:
            return x, y

    def mark_track(self, ax, xlims, ylims, masslbl=True, vvclbl=False):
        # Code responsible for making marks on the track and coloring it:
        #--------------------------------------------------------------------------------------
        # for placement on first loop if masking SGB+ (phase label 2):
        x = self.x
        y = self.y
        phasemask = self.phasemask
 
        if 2 in phasemask and False:
            past_growth1 = False
            for idx in range(len(x)):
                if idx > 10:
                    if y[idx] < y[idx - 1]:
                        # starting 1st decay:
                        past_growth1 = True
                    else:
                        continue

                if past_growth1:
                    # if growing again:
                    if y[idx] > y[idx - 1]:
                        loop1_idx = idx
                        break
                    else:
                        continue
            ax.scatter(x[idx], y[idx], marker='+')
            txtx = x[loop1_idx]
            txty = y[loop1_idx]
            xfac = 0.01
            yfac = 0.05

        # or else label at the bluest point.
        else:
            txtx = x[np.where(x == x.max())]
            xfac = 0.01
            txty = y[np.where(x == x.max())]
            yfac = 0.05

        # Label v/vcrit of the track:
        if vvclbl:
            txty = txty-abs(txty)*yfac
            ax.scatter(txtx, txty, marker = 'x', s=4, label=r'$\frac{\Omega}{\Omega_{crit}} = $' + "{:.1f}".format(self.rot),
                        c = [0.4+self.rot, 0, 0.6-self.rot], lw= 0.8, zorder = 9999)  

        # Label the mass of the track:
        #if self.minit not in plotted_masses:
        if masslbl:
            txtx = txtx+abs(txtx)*xfac
            txty = txty-abs(txty)*yfac
            if not (self.minit).is_integer():
                ax.text(txtx, txty, "{:.2f} ".format(self.minit) + r'$M_{\odot}$', fontsize = 16)
            else:
                ax.text(txtx, txty, "{:d} ".format(int(self.minit)) + r'$M_{\odot}$', fontsize = 16)

        # expand the axes to make space for new text.
        xlims, ylims = expand_lims(ax, x=txtx, y=txty, reverse_x=True)

        return xlims, ylims

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

def color_lineby_z(x, y, z, zmin, zmax, ax, cmap=plt.get_cmap('cool'), **kwargs):
        
    """
        See http://matplotlib.org/examples/pylab_examples/multicolored_line.html as the source.

        Colors a matplotlib line (w/ given x, y points) by a third array, z.       

    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, cmap=cmap, norm=plt.Normalize(zmin, zmax))
    lc.set_array(z)
    #lc.set_clim(vmin=0.0,vmax=1.0)

    # Try setting some key word arguments for line:

    # alpha
    try:
        lc.set_alpha(kwargs['alpha'])
    except KeyError:
        pass

    # label
    try:
        lc.set_label(kwargs['label'])
    except KeyError:
        pass

    # line style
    try:
        lc.set_linestyle(kwargs['ls'])
    except KeyError:
        pass

    # add colored line to axis.
    ax.add_collection(lc)

    return ax, lc

def expand_lims(ax, x=None, y=None, xlims=None, ylims=None, reverse_x=False, reverse_y=False, xpad=0, ypad=0):

    """
        Expands the x and y limits of the given axis; e.g. if plotting on an axis with limits
    that are too small to encompass the given x or y data points, expand the maxmimum and minimum values.

        Always give an axis to operate on.
 
        Give x and/or y data to have this method consider expanding the current limits to incorporate any
    data that may lie outside of the current axis limits. If x or y limits are given this will NOT happen,
    the method will assume those limits are meant to be applied.

        So to expand axis to incorporate given data values:

        >> expand_lims(ax, xdata, ydata)

    (optional to include both x and y data, could be one or the other also).

        Additionally, may include padding values to extend the axes further, by some percentage amount. So
    doing:

        >> expand_lims(ax, xpad=0.01)

    would extend the x-axis by 1% of its current limit values (in both directions, neg. and pos.). Right now
    padding operates exclusively from expanding, so if pads are given, padding will occur without expansion;
    perhaps should change this.

    """

    # if x or y padding values (given as percentages) are greater than 0...
    if xpad > 0 or ypad > 0:
        
        # create list of pads to add to either end of an axis:
        xpad = mkpad(ax.get_xlim(), xpad)
        ypad = mkpad(ax.get_ylim(), ypad)

        # if no limits are given, get the current axis limits.
        if ylims == None:
            ylims = ax.get_ylim()
        if xlims == None:
            xlims = ax.get_xlim()

        # add the created pads to the axis limits, extending them.
        ylims = [limit + ypad[i] for i, limit in enumerate(ylims)]
        xlims = [limit + xpad[i] for i, limit in enumerate(xlims)]

        # set the new padded limits.
        ax.set_ylim(ylims)
        ax.set_xlim(xlims)

        # return new limits:
        return xlims, ylims

    # if plot limits aren't given (if they are, use them instead):
    if ylims == None:

        # if no data is given, set the limits as current axis limits:
        if isinstance(y, type(None)):
            ylims = ax.get_ylim()

        # or else if data is given...
        else:
            # get the largest of the current limits & given data
            # (e.g. b/c want to expand axis lims if there is data
            # outside of them):
            ymax = max(max(y), max(ax.get_ylim()))
            # do the same check for the data min:         
            ymin = min(min(y), min(ax.get_ylim()))
            # designates new limits and reverses axes if told to:
            ylims = [ymin, ymax] if not reverse_y else [ymax, ymin]

    # same as above, but for x axis:
    if xlims == None:
        if isinstance(x, type(None)):
            xlims = ax.get_xlim()
        else:
            xmax = max(max(x), max(ax.get_xlim()))
            xmin = min(min(x), min(ax.get_xlim()))
            xlims = [xmin, xmax] if not reverse_x else [xmax, xmin]

    # Set limits, expanding to the extent that encompasses the full dataset, 
    # or else staying at the current limits if they already encompass all data.
    ax.set_ylim(ylims)
    ax.set_xlim(xlims)

    # return new limits:
    return xlims, ylims

def mkpad(lims, percent):

    """
     pads are some percentage of the given (axis) limits.
    
     -- 0% would add no pad
     -- 50% (0.5) would create a pad that is 50% of 
        the given limit value.
    
     pads are meant to be added to the given limit
     after their creation here.
    """

    # creates a list of pads (one element for either end of axis).
    pads = []
    for lim in lims:
        # if on the larger end of axis, create pads that will further
        # extend this end of the axis when added:
        if lim == max(lims):
            pads.append(abs(lim)*percent)
        # or else if on the other end of the axis, create pads that will
        # extend towards the other (negative) direction.
        elif lim == min(lims):
            pads.append(-abs(lim)*percent)
 
    return pads

def convert_name(aname):

    # For aesthetics: convert a given name (string) into a (matplotlib) plot friendly name.

    # If code is using TIGS photometry name convetions, then using 'TMASS'
    # not 2MASS, so...
    if 'TMASS' in aname:
        aname = aname.replace('TMASS', '2MASS')
    if aname == '2MASS_K':
        aname = '2MASS_Ks'

    # Keys = input names and values = the conversions:
    name_dictionary = {'2MASS_J': 'J',
                       '2MASS_Ks': 'K_s',
                       'Tycho_B': 'B_T',
                       'Tycho_V': 'V_T',
                       'Bessell_B': 'B',
                       'Bessell_V': 'V',
                       'log_Teff': 'log T_{eff}',
                       'log_L': 'log L/L_{\odot}'
                       }

    try:
        return name_dictionary[aname]

    except KeyError:
        # given name is not simply one of the keys...so it is trickier to convert or doesn't have a defined conversion.
        # condition to handle case where naem is a mathematical expression of valid names (e.g. 2MASS_J-2MASS_Ks)        
        for key in name_dictionary.keys():
            if key in aname:
                # replace the occurence of the name that should be converted:
                aname = aname.replace(key, name_dictionary[key])

        # If the given name ultimately doesn't contain any keys from the above dictioanry, the given name will be returned.

        return aname
