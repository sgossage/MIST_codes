import glob
import os

def isopath(vvcrit=None, cov=None, exttag=None):

    """
        This will find .iso files in the path structure defined below (relevant to my system).
    The inputs are a metallicity and a vvcrit value, but they are optional.
    """

    # Directory leading to the stored .iso files (for MY environment):
    storedir = os.path.join(os.environ['STORE_DIR'], 'MIST_v1.0', 'output')
   
    # Create the string used to search for desired iso files:

    search_str = ''

    if vvcrit is not None:
        search_str = search_str + '*_vvcrit{:.1f}'.format(vvcrit)
    if cov is not None:
        if cov != 0.016:
            search_str = search_str + '_cov{:.3f}'.format(cov)
        else:
            pass
    if exttag is not None:
        search_str += '_{:s}'.format(exttag)

    search_str = search_str + '*full.iso'

    if cov:
        if exttag is not None:
            pathtoisos = os.path.join(storedir, '*_vvcrit{:.1f}_cov{.3f}_{:s}/isochrones/{:s}'.format(vvcrit, cov, exttag, search_str))
        else:
            pathtoisos = os.path.join(storedir, '*_vvcrit{:.1f}_cov{.3f}/isochrones/{:s}'.format(vvcrit, cov, search_str))
    else:
        if exttag is not None:
            pathtoisos = os.path.join(storedir, '*_vvcrit{:.1f}_{:s}/isochrones/{:s}'.format(vvcrit, exttag, search_str))
        else:
            pathtoisos = os.path.join(storedir, '*_vvcrit{:.1f}/isochrones/{:s}'.format(vvcrit, search_str))
    #print(pathtoisos)
    # Returns a list of paths to each .iso file found in path structure above:
    isofilepaths = glob.glob(pathtoisos)
    #print(isofilepaths)
    
    return isofilepaths
