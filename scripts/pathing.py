import glob
import os

def isopath(vvcrit=None, cov=None):

    """
        This will find .iso files in the path structure defined below (relevant to my system).
    The inputs are a metallicity and a vvcrit value, but they are optional.
    """

    # Directory leading to the stored .iso files (for MY environment):
    storedir = os.environ['STORE_DIR']
   
    # Create the string used to search for desired iso files:

    search_str = ''

    if vvcrit:
        search_str = search_str + '*_vvcrit{:.1f}'.format(vvcrit)
    if cov:
        if cov != 0.016:
            search_str = search_str + '_cov{:.3f}'.format(cov)
        else:
            pass

    search_str = search_str + '*full.iso'

    if cov:
        pathtoisos = os.path.join(storedir, '*_vvcrit{:.1f}_cov{.3f}/isochrones/{:s}'.format(vvcrit, cov, search_str))
    else:
        pathtoisos = os.path.join(storedir, '*_vvcrit{:.1f}/isochrones/{:s}'.format(vvcrit, search_str))
    #print(pathtoisos)
    # Returns a list of paths to each .iso file found in path structure above:
    isofilepaths = glob.glob(pathtoisos)
    #print(isofilepaths)
    
    return isofilepaths
