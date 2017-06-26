# Aaron Dotter's grav. darkening code:
from pylab import *
from scipy.interpolate import RectBivariateSpline
from copy import copy
import os

#this creates two 2-D interpolants for C_T and C_L
def create_interpolants(npz):
    data=load(npz)
    C_L=data['C_L']
    C_T=data['C_T']
    omega=data['omega']
    inclination=data['inclination']
    f_T=RectBivariateSpline(y=omega, x=inclination, z=C_T)
    f_L=RectBivariateSpline(y=omega, x=inclination, z=C_L)
    return f_T, f_L

# angle in radians [0,pi/2]
def gdiso(iso, angle):
    #read in an isochrone set
    #fn=#'/n/conroyfs1/sgossage/MIST_grids/MIST_v1.0/output/feh_p0.15_afe_p0.0_vvcrit0.5/isochrones/MIST_v1.0_feh_p0.15_afe_p0.0_vvcrit0.5_full.iso'
    #x=ISO(filename) #ISO('MIST_v1.0_feh_p0.00_afe_p0.0_vvcrit0.4_full.iso')

    f_T, f_L = create_interpolants(os.path.join(os.environ['MIST_CODE_DIR'], 'GD.npz'))

    #choose one isochrone for demonstration
    #iso=x.isos[x.age_index(8.5)]
    i=squeeze(where( (iso['EEP']<9999) & (iso['EEP']>0) ))#<808 & > 202
    origiso = copy(iso)
    # all cols where eep reqs specified above are met:
    iso=iso[:][i]
    omega=iso['surf_avg_omega_div_omega_crit']
    T=pow(10,iso['log_Teff'])
    L=pow(10,iso['log_L'])

    #figure(1)
    #plot the intrinsic model quantities
    #plot(T,log10(L), color='Lime', label='Intrinsic')

    #plot the gravity-darkened models at i=0, "edge-on"
    if angle == 0.0: 
        C_T=f_T.ev(yi=omega, xi=zeros(len(omega)))
        C_L=f_L.ev(yi=omega, xi=zeros(len(omega)))
    #plot(C_T*T, log10(C_L*L), color='Blue', label=r'$i=0$')

    #plot the gravity-darkened models at i=90 deg, "face-on"
    else:
        C_T=f_T.ev(yi=omega, xi=(angle)*ones(len(omega)))
        C_L=f_L.ev(yi=omega, xi=(angle)*ones(len(omega)))

    gdlTeff = log10(C_T*T)
    gdlLum = log10(C_L*L)

    #print(gdlTeff == log10(T))

    #tempT = copy(origiso['log_Teff'])
    for j, idx in enumerate(i):
        origiso['log_Teff'][idx] = gdlTeff[j]
        origiso['log_L'][idx] = gdlLum[j]
        #print(gdlLum[j])

    #print(tempT == origiso['log_Teff'])

    return origiso

# angle in radians [0,pi/2]
def gdtrack(track, angle):
    #read in an isochrone set
    #fn=#'/n/conroyfs1/sgossage/MIST_grids/MIST_v1.0/output/feh_p0.15_afe_p0.0_vvcrit0.5/isochrones/MIST_v1.0_feh_p0.15_afe_p0.0_vvcrit0.5_full.iso'
    #x=ISO(filename) #ISO('MIST_v1.0_feh_p0.00_afe_p0.0_vvcrit0.4_full.iso')

    f_T, f_L = create_interpolants(os.path.join(os.environ['MIST_CODE_DIR'], 'GD.npz'))

    #choose one isochrone for demonstration
    #iso=x.isos[x.age_index(8.5)]
    i=squeeze(where( (track['phase']<9999) & (track['phase']>-1) ))#<808 & > 202
    origtrack = copy(track)
    # all cols where eep reqs specified above are met:
    track=track[:][i]
    omega=track['surf_avg_omega_div_omega_crit']
    T=pow(10,track['log_Teff'])
    L=pow(10,track['log_L'])

    #figure(1)
    #plot the intrinsic model quantities
    #plot(T,log10(L), color='Lime', label='Intrinsic')

    #plot the gravity-darkened models at i=0, "edge-on"
    if angle == 0.0: 
        C_T=f_T.ev(yi=omega, xi=zeros(len(omega)))
        C_L=f_L.ev(yi=omega, xi=zeros(len(omega)))
    #plot(C_T*T, log10(C_L*L), color='Blue', label=r'$i=0$')

    #plot the gravity-darkened models at i=90 deg, "face-on"
    else:
        C_T=f_T.ev(yi=omega, xi=(angle)*ones(len(omega)))
        C_L=f_L.ev(yi=omega, xi=(angle)*ones(len(omega)))

    gdlTeff = log10(C_T*T)
    gdlLum = log10(C_L*L)

    #print(gdlTeff == log10(T))

    #tempT = copy(origiso['log_Teff'])
    for j, idx in enumerate(i):
        origtrack['log_Teff'][idx] = gdlTeff[j]
        origtrack['log_L'][idx] = gdlLum[j]
        #print(gdlLum[j])

    #print(tempT == origiso['log_Teff'])

    return origtrack

def rw_gdiso(filename, angle):
    # re-writing iso file with grav. darkening effects:
    
    # convert rad to deg
    angle_deg = (180./pi)*angle
    outfname = filename.split('.iso')[0] + '_gdark{:.1f}.iso'.format(angle_deg)
    # Check if the .iso file already exists:
    if os.path.isfile(outfname):
        print("The .iso file \"{:s}\" already exists and will be used.".format(outfname))
        return outfname
    # If it doesn't, create it:
    else:

        isoobj = ISO(filename)

        numisos = len(isoobj.isos)

        # organizing iso data blocks:
        isoblocks = []
        print('Organizing iso data...')
        for i in range(numisos):
            block = []
            numrows = len(isoobj.isos[i]['EEP'][:])
            iso = copy(gdiso(isoobj.isos[i], angle))

            for r in range(numrows):
                row = []
                for name in isoobj.hdr_list:
                    if name != 'EEP':
                        # adopting precision in MIST .iso files:
                        row.append("{:.16E}".format(iso[name][r]))
                    else:
                        # EEPs are recorded as ints:
                        row.append("{:d}".format(iso[name][r]))

                block.append(row)

            isoblocks.append(block)
        
        # write data to .iso file in format similar to A. Dotter's iso code.        
        print("Writing {:s}...".format(outfname))
        with open(outfname, 'w+') as outf:
            outf.write('# MIST version number  = {:s}\n'.format(isoobj.version['MIST']))
            outf.write('# MIST revision number = {:s}\n'.format(isoobj.version['MESA']))
            outf.write('# --------------------------------------------------------------------------------------\n')
            outf.write('#{:>7s}{:>14s}{:>9s}{:>9s}'.format(isoobj.abun.keys()[1],isoobj.abun.keys()[0],isoobj.abun.keys()[2],isoobj.abun.keys()[3]))
            outf.write('{:>9s}\n'.format('v/vcrit'))
            outf.write('#{:>7s}'.format("{:.4f}".format(isoobj.abun.values()[1]))) # Yinit
            outf.write('{:>14s}'.format(eformat(isoobj.abun.values()[0], 5, 2))) # Zinit
            outf.write('{:>9s}'.format("{:.2f}".format(isoobj.abun.values()[2]))) # [Fe/H]
            outf.write('{:>9s}'.format("{:.2f}".format(isoobj.abun.values()[3]))) # [a/Fe]
            outf.write('{:>9s}\n'.format("{:.2f}".format(isoobj.rot))) #v/vcrit
            outf.write('# --------------------------------------------------------------------------------------\n')
            outf.write('# number of isochrones =   {:d}\n'.format(len(isoblocks)))
            outf.write('# --------------------------------------------------------------------------------------\n')
            # write blocks of iso data:
            for i, isoblock in enumerate(isoblocks):
                # for a block, write num EEPs, num columns
                outf.write('# number of EEPs, cols =   {:d}   {:d}\n'.format(len(isoblock), len(isoblock[0])))
                outf.write('#{:>4s}'.format('1'))
                for i in range(2,80):
                    outf.write('{:>32s}'.format(str(i)))
                outf.write('\n')

                # for a block, write header:
                for name in isoobj.hdr_list:
                    if name == 'EEP':
                        outf.write('#{:>4s}'.format(name))
                    else:
                        outf.write('{:>32s}'.format(name))

                # for an iso block, write rows
                for blockrow in isoblock:
                    outf.write('\n')
                    for j, val in enumerate(blockrow):
                        if j == 0:
                            outf.write(' {:>4s}'.format(str(val)))
                        else:
                            outf.write('{:>32s}'.format(eformat(float(val), 16, 3)))

                if i < len(isoblocks) - 2:
                    outf.write('\n\n\n')

        return outfname

def rw_gdeep(eep, filename, angle):
    # re-writing iso file with grav. darkening effects:
    
    # convert rad to deg
    angle_deg = (180./pi)*angle
    outfname = filename.split('.track.eep')[0] + '_gdark{:.1f}.track.eep'.format(angle_deg)
    # Check if the .iso file already exists:
    if os.path.isfile(outfname):
        print("The .iso file \"{:s}\" already exists and will be used.".format(outfname))
        return outfname
    # If it doesn't, create it:
    else:

        eepobj = eep #EEP(filename)

        numeeps = len(eepobj.eeps)

        # organizing iso data blocks:
        trackblock = []
        print('Organizing track data...')
        numrows = numeeps
        track = copy(gdtrack(eepobj.eeps, angle))

        for r in range(numrows):
            row = []
            for name in eepobj.hdr_list:
                # adopting precision in MIST .eep files:
                row.append("{:.16E}".format(track[name][r]))

            trackblock.append(row)

        # get the content preceding data...
        with open(filename) as origf:
            origcont = origf.readlines()[0:10]
        
        # write data to .iso file in format similar to A. Dotter's iso code.        
        print("Writing {:s}...".format(outfname))
        with open(outfname, 'w+') as outf:
            for line in origcont:
                outf.write(line)

            # write track data:
            # for a block, write col number:
            outf.write('#{:>4s}'.format('1'))
            for i in range(2,80):
                outf.write('{:>32s}'.format(str(i)))
            outf.write('\n')

            # write header:
            for name in eepobj.hdr_list:
                if name == eepobj.hdr_list[0]:
                    outf.write('#{:>4s}'.format(name))
                else:
                    outf.write('{:>32s}'.format(name))

            # write data rows
            for blockrow in trackblock:
                outf.write('\n')
                for j, val in enumerate(blockrow):
                    if j == 0:
                        outf.write(' {:>4s}'.format(str(val)))
                    else:
                        outf.write('{:>32s}'.format(eformat(float(val), 16, 3)))

            outf.write('\n\n\n')

        return outfname

def eformat(f, prec, exp_digits):
    s = "%.*e"%(prec, f)
    mantissa, exp = s.split('e')
    # add 1 to digits as 1 is taken by sign +/-
    return "%sE%+0*d"%(mantissa, exp_digits+1, int(exp))

#plot(C_T*T, log10(C_L*L), color='Red', label=r'$i=\pi/2$')

#xlim(11000,2000)
#ylim(0,3.0)

#xlabel(r'$T_{eff} (K)$', fontsize=18)
#ylabel(r'$L/L_{\odot}$', fontsize=18)

#legend(loc='lower right', fontsize=18, fancybox=True)

#show()
