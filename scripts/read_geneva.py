import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

class star:

    """
        This class holds parameters that correspond to a Geneva (cite) stellar model; i.e. a model of a particular mass.

        Useage:
        -------
            Invoke with
 
                star(filename), 

            where filename pertains to a Geneva stellar evolution track file downloaded from the Geneva project site. Or 
            else, 

                star(vvc, minit)
 
            where vvc is the v/vcrit (rotation) value of the model, and minit is the desired mass for the model. Use
            e.g.
          
                star.hdr_list

            to print a list of the available data columns extracted from the model file. Then, either

                star.data[column_name]

             or
 
                star.dataframe[column_name]

             to retrieve a particular data column from the file.

        Arguments:
        ----------
            Either supply a filename through fname, or else a vvc and minit; in the latter case, the model will be sought
            out in a directory designated below by the variable 'fpath'.

            + fname (str): Path to the model file
            + vvc (float): v/vcrit value of the model; supply it as a float between 0 and 1.0
            + minit (float): Mass of the model; also a float; will be truncated to 2 digit precision.
    """

    def __init__(self, fname=None, vvc=None, minit=None):
        # If given a filename, use it:
        if fname is not None:
            # extract the v/vcrit and initial model mass from the filename:
            self.fname = fname
            self.minit = fname.split('Z')[0].split('M')[-1]
            if 'p' in self.minit:
                self.minit = float('.'.join(self.minit.split('p')))
            else:
                self.minit = float(self.minit)

            self.vvc = (self.fname).split('V')[-1].split('.dat')[0]
            if self.vvc == 'c':
                self.vvc = 1.0
            else:
                self.vvc = float('{:.1f}'.format(float(self.vvc)/10.0))

        # Or else without a filename given, use the supplied v/vcrit and initial mass values
        # to find the desired file:
        else:
            # Form the initial mass string corresponding to Geneva file name conventions:
            minitstr = '{:.2f}'.format(minit)
            minitstr = minitstr.split('.')
            # Geneva filenames use MXXX to distinguish the model's mass, and MXpXX for decimal numbers.
            # Decimals:
            if float(minitstr[-1]) > 0.0:
                # if there is a non-zero value following decimal, join via 'p'; accordance with the Geneva file name conventions.
                minitstr = 'p'.join((minitstr[0], minitstr[1]))
            # Integers:
            else:
                s = ['0','0','0']
                digits = minitstr[0].split()
                for i, digit in enumerate(digits):
                    s[-1+i] = digits[-1+i]

                minitstr = "".join(s)

            # initial mass of the star:
            self.minit = minit
            # v/vcrit of star (given by user) must be btwn 0 and 1:
            assert 0.0 <= vvc <= 1.0, "Critical rotation value must be between 0 and 1."
            self.vvc = vvc

            # Geneva filenames have V/Vcrit string corresp. to VX; X = [0,9] or 'c'  when X=1 -- demarcating v/vc = 0.0, 0.1, 0.6, 1.0, etc.
            vvcstr = int(vvc*10)
            # critical rotation:
            if vvcstr == 10:
                vvcstr = 'c'
            # sub-critical:
            else:
                vvcstr = str(vvcstr)
                
            # my path to Geneva files (**system dependent**):
            fpath = '/home/seth/Research/Geneva'

            # Resultant file name:
            self.fname = os.path.join(fpath, 'z0.014','tracks','dense','M{:s}Z14V{:s}.dat'.format(minitstr, vvcstr))

        # Call read_file() to extract model data table and header labels:
        self.hdr_list, self.hdr_units, self.data = self.read_file()
        # Organize data table into a pandas data frame:
        self.df = pd.DataFrame(self.data, columns=self.hdr_list)

        # Label for plots:
        self.lbl='Geneva: M = {:.2f}'.format(self.minit)+r' $M_{\odot}$, $\Omega/\Omega_c$'+' = {:.1f}'.format(self.vvc)

    def read_file(self):
        
        """

            Reads in the Geneva track file. Arguments are not needed as this uses the class' filename property
        which was created during the initialization process.

        Returns:
        --------
            This returns the data column header names, their units, and a dictionary with column names as keys and
        corresponding data columns as items.
                
        """
        
        # Extract data columns, skip header lines:
        data = np.genfromtxt(self.fname, skip_header=2)

        # Read the header lines:
        print("Reading in {:s}".format(self.fname.split('/')[-1]))
        with open(self.fname) as f:
            hdr_lines = f.readlines()[:2]
            # Labels/names:
            hdr_list = hdr_lines[0].split()
            # Units:
            hdr_units = hdr_lines[1].split()

        # Dictionary to organize data columns as items w/ column labels as keys:
        datadict = dict(zip(hdr_list, data.T))
        
        # Return header info & data dictionary:
        return hdr_list, hdr_units, datadict

    def plot_HR(self, ax, shade=None, legloc=None, label=False, masslbl=False, vvclbl=False, show=False, **kwargs):

        # Get the log L and log Teff for an HR diagram.    
        geny = logL = self.data['lg(L)']
        genx = logTeff = self.data['lg(Teff)']

        self.x = genx
        self.y = geny

        # Create axes if necc:
        if ax == None:
            fig = plt.figure(figsize=(16,9))
            ax = fig.add_subplot(111)

        if shade != None:
            lc = plt.cm.Dark2(shade)
            kwargs['c'] = lc

        if label:
            kwargs['label'] = self.lbl

        base_line, = ax.plot(genx, geny, **kwargs)

        if masslbl or vvclbl:
            self.mark_track(ax, masslbl=masslbl, vvclbl=vvclbl)

        if legloc != None:
            ax.legend(loc = legloc, prop={'size':8})

        if show:
            plt.show()

        return genx, geny, base_line.get_color() 

    def mark_track(self, ax, masslbl=True, vvclbl=False):
        # Code responsible for making marks on the track and coloring it:
        #--------------------------------------------------------------------------------------
        # for placement on first loop if masking SGB+ (phase label 2):
        x = self.x
        y = self.y
        txtx = x[np.where(x == x.max())]
        xfac = 1.003
        txty = y[np.where(x == x.max())]
        yfac = 0.997

        # Label v/vcrit of the track:
        if vvclbl:
            ax.scatter(txtx, txty*yfac, marker = 'x', s=4, label=r'$\frac{\Omega}{\Omega_{crit}} = $' + "{:.1f}".format(self.rot),
                        c = [0.4+self.rot, 0, 0.6-self.rot], lw= 0.8, zorder = 9999)  

        # Label the mass of the track:
        #if self.minit not in plotted_masses:
        if masslbl:
            ax.text(txtx*xfac, txty*yfac, "{:.2f} ".format(self.minit) + r'$M_{\odot}$', fontsize = 10)

        return

"""
def plt_abun(vvc):

    """"""
       Plots the ratio of intial surface abundances of C12, N14, O16 compared to their
    end of MS values for Geneva and MIST stars @ a particular v/vcrit value.

        + Rotational mixing is expected to increase the amount of n14
        + Decrease the amount of C12
        + decrease the amount of O16

    This is meant to show how MIST rotational mixing compares to Geneva's. Specialized and maybe unnecessary.

    """"""
    assert 0.0 <= vvc <= 1.0, "Given v/vcrit must be between 0 and 1."

    surfc12_ratio = []
    surfo16_ratio = []
    surfn14_ratio = []

    mist_surfc12_ratio = []
    mist_surfo16_ratio = []
    mist_surfn14_ratio = []

    imasses = []

    allbutMS_mask = [-1, 2, 3, 4, 5]
    MIST_massgrid = np.concatenate((np.arange(1,2,0.04),np.arange(2,8.2,0.2)))
    #print(MIST_massgrid)

    for i, mass in enumerate(MIST_massgrid):

        # get a MIST model (Z=0.014 [Fe/H]=0.0 is all I have for Geneva atm:
        miststar = rmm.EEP(feh = 0.00, vvcrit = vvc, mass = mass)

        # try to get a corresponding Geneva star:
        try:
            genevastar = star(vvc=vvc, minit=mass)
            surfc12_ratio.append(np.log10(genevastar.df['12C_surf'].values[109] / genevastar.df['12C_surf'].values[1]))
            surfo16_ratio.append(np.log10(genevastar.df['16O_surf'].values[109] / genevastar.df['16O_surf'].values[1]))
            surfn14_ratio.append(np.log10(genevastar.df['14N_surf'].values[109] / genevastar.df['14N_surf'].values[1]))
            # a list of all surface abundances given below during the MS phase of the MIST star:
            mist_surfabun = miststar.get_masked(['surface_c12', 'surface_o16', 'surface_n14'], phasemask = allbutMS_mask)
            mist_surfc12_ratio.append(np.log10(mist_surfabun[0][-1] / mist_surfabun[0][0]))
            mist_surfo16_ratio.append(np.log10(mist_surfabun[1][-1] / mist_surfabun[1][0]))
            mist_surfn14_ratio.append(np.log10(mist_surfabun[2][-1] / mist_surfabun[2][0]))
            # the initial masses
            imasses.append(mass)

        # in case the Geneva star is not recoverable:
        except IOError:
            print('Geneva star: v/vcrit = {:.1f}, Mass = {:.2f}'.format(vvc, mass))


    geneva_surfabun = {'surfc12': surfc12_ratio, 'surfo16': surfo16_ratio, 'surfn14': surfn14_ratio}
    mist_surfabun = {'surfc12': mist_surfc12_ratio, 'surfo16': mist_surfo16_ratio, 'surfn14': mist_surfn14_ratio}

    # returns the above dictionaries; use them to create plots.
    # E.g. ---
    # plt.plot(imasses, geneva_surfabun['surfc12'])
    # plt.plot(imasses, mist_surfabun['surfc12'])
    # ...

    return geneva_surfabun, mist_surfabun, imasses
"""

    

        


