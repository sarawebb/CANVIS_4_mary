"""location_canvis.py -- Input an RA and DEC cut out postage stamps around associated RA and DEC for all available data for given field; Input also the DWF_run, which specifies where CANVIS outputs will be saved. 

Usage: location_canvis [-h] [-v] [--debug] <RA> <DEC> <field> <DWF_run>

Arguments:
    RA (float)
        The RA in degrees.
    DEC (float)
        The DEC in degrees.
    field (string)
        The DWF field name. 
    DWF_run (string)
        The DWF run date/ name. This specifies the folder under which gifs are saved at: /fred/oz100/CANVIS/cand_images/DWF_run

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]     
    --debug                                 Output more for debugging [default: False]

Example:
    python location_canvis.py -v 121.0000000 -78.2600000 8hr zhangtesttest
    python location_canvis.py -v 125.0000000 -78.2600000 8hr zhangtesttest
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


#IDs = np.arange(1,5)
#field = '8hr'
#run ='SARATEST'


###################################### functions #################################



def RAdec_to_RAsex(fRAdec):
    fratotsec = (math.fabs(float(fRAdec))*3600.0)
    frah2 = (math.modf(fratotsec/3600.0)[1])
    fram2 = (math.modf((fratotsec-(frah2*3600.0))/60.0)[1])
    fras2 = (fratotsec-(frah2*3600.0)-(fram2*60.0))
    if round(fras2, 2) == 60.00:
        fram2 = fram2 + 1
        fras2 = 0
        if round(fram2, 2) == 60.00:
            frah2 = frah2 + 1
            fram2 = 0
    if round(fram2, 2) == 60.00:
        frah2 = frah2 + 1
        fram2 = 0
    if int(frah2) == 24 and (int(fram2) != 0 or int(fras2) != 0):
        frah2 = frah2 - 24
    fRAsex = '%02i' % frah2 + ' ' + '%02i' % fram2 + ' ' + ('%.3f' % float(fras2)).zfill(6)
    return fRAsex


def DEdec_to_DEsex(fDEdec):
    fdetotsec = (math.fabs(float(fDEdec))*3600.0)
    fded2 = (math.modf(fdetotsec/3600.0)[1])
    fdem2 = (math.modf((fdetotsec-(fded2*3600.0))/60.0)[1])
    fdes2 = (fdetotsec-(fded2*3600.0)-(fdem2*60.0))
    if float(fDEdec) < 0:
        fded2sign = '-'
    else:
        fded2sign = '+'
    fDEsex = fded2sign + '%02i' % fded2 + ' ' + '%02i' % fdem2 + ' ' + ('%.2f' % float(fdes2)).zfill(5)
    return fDEsex


def RAsex_to_RAdec(fRAsex):
    frah = float(fRAsex[0:2])
    #print(frah)
    fram = float(fRAsex[2:4])
    #print(fram)
    fras = float(fRAsex[4:])
    #print(fras)
    #fRAdec = (frah*3600.0+fram*60.0+fras)/3600.0
    fRAdec = ((1/1 * frah) + (1/60 *fram) + (1/3600 *fras))* (360/24)
    return fRAdec
    
      
def DEsex_to_DEdec(fDEsex):
    fded = float(fDEsex[0:3])
    print(fded)
    fdem = float(fDEsex[3:5])
    print(fdem)
    fdes = float(fDEsex[5:])    
    print(fdes)
    fDEdec = (math.fabs(fded)*3600.0+fdem*60.0+fdes)/3600.0
    if fDEsex[0] == '-':
        fDEdec = fDEdec * -1
    return fDEdec


def location_canvis(RA, DEC, field,run,verbose=False,debugmode=False):
    print('\n#########################')
    print('#  CANVIS HAS STARTED   #')
    print('#########################\n')
    print(f'CANVIS will make cutout for location RA: {RA} DEC: {DEC}')

    '''RA, DEC, CANVIS will go through all the data on this field  
    and extract postage stamps around the given RA and DEC.'''
    print("CANVIS will extract postage stamps around RA %s DEC %s for field %s." %(RA, DEC,field))
    path ='/fred/oz100/pipes/DWF_PIPE/MARY_WORK/'+field+'_*_*_*/*/images_resampled/sci_*.resamp.fits'
    fitsfileslist = glob.glob(path)
    mydic = SortedDict()
    for i  in fitsfileslist:
        with fits.open(i) as hdu:
            size = 200 
            w = WCS(hdu[0].header)
            head = hdu[0].header
            date = dt.datetime.strptime(head['DATE'], '%Y-%m-%dT%H:%M:%S')
            xlim=head['NAXIS1']
            ylim=head['NAXIS2']
            pixcrd_im = np.array([[xlim, ylim]], np.float_)
            world_im = w.wcs_pix2world(pixcrd_im, 1)
            pixx_im, pixy_im = world_im[0][0], world_im[0][1]
            corners=w.calc_footprint()

            corner_1 = corners[0]
            corner_2 = corners[1]
            corner_3 = corners[2]
            corner_4 = corners[3]
            differnce = corner_1 - corner_2 
                
            pixcrd = np.array([[RA, DEC]], np.float_)
            worldpix = w.wcs_world2pix(pixcrd, 1)
            pixx, pixy = worldpix[0][0], worldpix[0][1]

            if float( corner_4[0]) <= float(RA) <=float(corner_1[0]) and float(corner_2[1]) >= float(DEC) >= float(corner_1[1]):
                path = i
                mydic[date] =[path, pixx, pixy]
                if debugmode:
                    print(mydic)

    images_found=False
    for i, (key, (path, pixx, pixy)) in enumerate(mydic.items()):
        images_found=True
        if debugmode == True:
            print('variable run is: {} and the type is: {}'.format(run,type(run)))
            print('variable RA is: {} and the type is: {}'.format(RA,type(RA)))
            print('variable DEC is: {} and the type is: {}'.format(DEC,type(DEC)))
            print('variable field is: {} and the type is: {}'.format(field,type(field)))
        path_cand = '/fred/oz100/canvis/cand_images/'+ run + '/cand_'+RA+DEC+'_'+ field +'_'+ run +'/'
        path_cutout = '/fred/oz100/CANVIS/cand_images/'+ run +'/cand_'+RA+DEC+'_'+ field +'_'+ run +'/'+RA+DEC+'_'+run+'_cutout_'+format(i, '03')  
        if not os.path.exists(path_cand):
            os.makedirs(path_cand, 0o755)
        else:
            pass
        if not os.path.exists(path_cutout):
            os.makedirs(path_cutout, 0o755)
        else:
            pass
        size = 200
        with fits.open(path) as hdu:
            nom_data = (hdu[0].data - np.min(hdu[0].data))/(np.max(hdu[0].data)-np.min(hdu[0].data))
            cutout = Cutout2D(hdu[0].data, (pixx, pixy), size, wcs= w)
            hdu[0].data = cutout.data
            hdu[0].header['CRPIX1'] = cutout.wcs.wcs.crpix[0]
            hdu[0].header['CRPIX2'] = cutout.wcs.wcs.crpix[1]
            hdu.writeto(path_cutout+'.fits', overwrite = True)
            plt.axis('off')
            plt.imshow(hdu[0].data, cmap='gray')
            plt.colorbar()
            plt.savefig(path_cutout+'.png', overwite=True)
            plt.close()

    if images_found:
        files = []	
        path_cutout = '/fred/oz100/CANVIS/cand_images/'+ run +'/cand_'+RA+DEC+'_'+ field +'_'+ run +'/'
        for cutouts in os.listdir(path_cutout):
            if cutouts.endswith('.png'):
                files.append(path_cutout + cutouts)
        writer = imageio.get_writer(str(path_cutout)  + '_VIDEO.gif', fps =3)
        video_loc = str(path_cutout)  + '_VIDEO.gif'
        for i in files:
            writer.append_data(imageio.imread(i))
        writer.close()
    else:
        print('\nCANVIS did not find any images to create a gif with, sorry!\n')

    if images_found:
        print('\nDone! Look for your outputs here: /fred/oz100/CANVIS/cand_images/%s' % run)

    print('\n###########################')
    print('#  CANVIS HAS Finished    #')
    print('#  Enjoy and discover!    #')
    print('###########################\n')

    return video_loc


if __name__ == "__main__":

    # Read in arguments
    arguments = docopt.docopt(__doc__, options_first=True)

    # Mandatory arguments
    RA = arguments['<RA>']
    DEC = arguments['<DEC>']
    field = arguments['<field>']
    run   = arguments['<DWF_run>']
    print(arguments)
    # Optional arguments
    verbose = arguments['--verbose']
    debugmode       = arguments['--debug']
    
    if debugmode:
        print(arguments)

    location_canvis(RA,DEC,field,run,verbose=verbose,debugmode=debugmode)
