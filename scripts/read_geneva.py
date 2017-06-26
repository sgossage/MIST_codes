import numpy as np
import pandas as pd

class star:

    def __init__(self, fname):
        self.fname = fname

        self.hdr_list, self.data = self.read_file()


    def read_file(self):
        
        """

        Reads in the EEP file.
        
        Args:
            filename: the name of .track.eep file.
                
        """
        
        data = np.genfromtxt(self.fname, names=True) #np.genfromtxt(self.fname, names=True)

        hdr_list = data.dtype.names
        
        return hdr_list, data
