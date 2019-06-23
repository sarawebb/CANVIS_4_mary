"""run_canvis.py -- Input a field, Mary ID and seed ID, cut out postage stamps around associated RA and DEC for all available data for given field and CCD; Input also the DWF_run, which specifies where CANVIS outputs will be saved. 

Usage: canvis_RA_persieve.py [-h] [-v] [--debug] <persieve_canvis_file> <run>

Arguments:
    persieve_canvis_file (string)
        The name of the persieve file that you want to run CANVIS on.
    run (string)
        Which DWF run is the data from. Eg. june19 

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]     
    --debug                                 Output more for debugging [default: False]

Example:
    python canvis_RT_persieve.py PerSieve_CANVIS_requests.data
 """

import docopt
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.wcs import WCS
from astropy.io import fits
import sys
import math
import os
import glob
import sys
from sortedcontainers import SortedDict
import datetime as dt
import imageio
import os
from PIL import Image
from matplotlib.colors import LogNorm
from astropy.nddata.utils import Cutout2D
from astropy import units as u
from PyAstronomy import pyasl
from run_canvis import run_canvis  

#### Import list of IDs requested by PerSieve ### 

def make_persieve_canvis_requests(persieve_canvis_file,run, verbose=False, debug=False):

    persieve_cat_path = '/fred/oz100/RT_CANVIS/PerSieve_Cat/PerSieve_CANVIS_requests.dat'
    IDs = np.loadtxt(persieve_cat_path, usecols= 0)
    fields = np.loadtxt(persieve_cat_path, usecols =1, dtype=str)
    print(IDs, fields)
    
    for ID, field in zip(IDs, fields):
        print(ID, field)
        run_canvis(field,ID, run)




if __name__ == "__main__":

    # Read in arguments
    arguments       = docopt.docopt(__doc__)

    # Mandatory arguments
    persieve_canvis_file = arguments['<persieve_canvis_file>']
    run = arguments['<run>']
    # Optional arguments
    verbose = arguments['--verbose']
   

    make_persieve_canvis_requests(persieve_canvis_file,run,verbose=verbose)
