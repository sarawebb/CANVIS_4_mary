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
from argparse import ArgumentParser
########################## Enter ID ################################## 

RA = 246.5000000   
DEC = -73.0000000
field = 'ngc6101' 
template = 'REQ16A_20160727_8c29d52-g-20160728T041310_si1_dither_seeing12_join.fits'

parser = ArgumentParser()

parser.add_argument('--RA', required=False, 
                        help='RA of the object you want to use CANVIS on.',
                        dest='RA', nargs='+')
			
parser.add_argument('--DEC', required=False, 
                        help='DEC of the object you want to use CANVIS on.',
                        dest='DEC', nargs='+')

parser.add_argument('--field', required=False, 
                        help='field of the object you want to use CANVIS on.',
                        dest='field', nargs='+')



###################################### funcations #################################



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


###################-------- CANdidate VISualisation  ---------#######################
path ='/fred/oz100/pipes/DWF_PIPE/TEMPLATES/' + field +'/'+ template

mydic = SortedDict()


with fits.open(path) as hdu:
    size = 1200
    w = WCS(hdu[0].header)
    #print(w)
    head = hdu[0].header
    date = dt.datetime.strptime(head['DATE'], '%Y-%m-%dT%H:%M:%S')
    xlim=head['NAXIS1']
    ylim=head['NAXIS2']
    print(xlim, ylim)
    pixcrd_im = np.array([[xlim, ylim]], np.float_)
    world_im = w.wcs_pix2world(pixcrd_im, 1)
    pixx_im, pixy_im = world_im[0][0], world_im[0][1]

    pixcrd = np.array([[RA, DEC]], np.float_)
    worldpix = w.wcs_world2pix(pixcrd, 1)
    pixx, pixy = worldpix[0][0], worldpix[0][1]
    
    cutout = Cutout2D(hdu[0].data, (pixx, pixy), size, wcs= w)
    hdu[0].data = cutout.data
    hdu[0].header['CRPIX1'] = cutout.wcs.wcs.crpix[0]
    hdu[0].header['CRPIX2'] = cutout.wcs.wcs.crpix[1]
    hdu.writeto('TEST.fits', overwrite = True)
    plt.axis('off')
    plt.imshow(hdu[0].data, cmap='gray', vmin= -1, vmax =10)
    plt.savefig('TEST.png', overwite=True)
    plt.close()
