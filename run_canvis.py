"""run_canvis.py -- Input a field, Mary ID and seed ID, cut out postage stamps around associated RA and DEC for all available data for given field and CCD; Input also the DWF_run, which specifies where CANVIS outputs will be saved. 

Usage: run_canvis [-h] [-v] [--debug] <field> <ID> <seed> <DWF_run>

Arguments:
    field (string)
        The DWF field name. 
    ID (integer)
        Mary ID of the object you want to use CANVIS on. Given a field, Mary IDs are unique.
    seed (string)
        The seed tells you PLEASEEDITSARA.
    DWF_run (string)
        The DWF run date/ name. This specifies the folder under which gifs are saved at: /fred/oz100/CANVIS/cand_images/DWF_run

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]     
    --debug                                 Output more for debugging [default: False]

Example:
    python ~/jlzhang/DWF_runaids/run_canvis.py -v 8hr 2 rt_TEST3 DWF_foo
    python run_canvis.py -v 8hr 2 rt_TEST3 DWF_foo
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


#ID = 2
#field = '8hr'
#seed ='rt_TEST3'
#run ='jlzhangtest'


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


def run_canvis(field,ID,seed,run,verbose=False,debugmode=False):
    print('\n#########################')
    print('#  CANVIS HAS STARTED   #')
    print('#########################\n')

    translist_path = '/home/fstars/MARY4OZ/transients/transients_coo_'+field+'_'+seed+'.txt'
    print('CANVIS will make a gif for field %s Mary ID number %s.\n' % (field,ID) )
    print('CANVIS will do this by reading %s\n' %translist_path)

    '''Look for the CCD number, RA and DEC of the Mary ID entry that matches inputs.'''
    with open(translist_path) as f:
        for line in f:
            line = line.split()
            #print(line[0])
            if int(line[0]) == ID:
                ra = str(line[1])
                dec = str(line[2])
                #field = str(line[3])           
                ccd_num = str(line[6])   

    if debugmode:
        print(ra,dec,ccd_num)

    '''Given the CCD number, RA, DEC, go through all the data on this field with this CCD 
    and extract postage stamps around the given RA and DEC.'''
    print("CANVIS will extract postage stamps around RA %s DEC %s for field %s on CCD %s for all seed IDs and dates." %(ra, dec,field,ccd_num))
    path ='/fred/oz100/pipes/DWF_PIPE/MARY_WORK/'+field+'_*_*_*/ccd' + ccd_num+'/images_resampled/sci_' + ccd_num+'.resamp.fits'
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
                
            pixcrd = np.array([[ra, dec]], np.float_)
            worldpix = w.wcs_world2pix(pixcrd, 1)
            pixx, pixy = worldpix[0][0], worldpix[0][1]

            if float( corner_4[0]) <= float(ra) <=float(corner_1[0]) and float(corner_2[1]) >= float(dec) >= float(corner_1[1]):
                path = i
                mydic[date] =[path, pixx, pixy]
                if debugmode:
                    print(mydic)

    for i, (key, (path, pixx, pixy)) in enumerate(mydic.items()):
        path_cand = '/fred/oz100/CANVIS/cand_images/'+ run + '/cand_'+format(ID, '05')+'_'+ field +'_'+ run +'/'
        path_cutout = '/fred/oz100/CANVIS/cand_images/'+ run +'/cand_'+format(ID, '05')+'_'+ field +'_'+ run +'/cand_'+format(ID, '05')+'_'+run+'_cutout_'+format(i, '03')  
        if not os.path.exists(path_cand):
            os.makedirs(path_cand, 0o755)
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

    files = []	
    path_cutout = '/fred/oz100/CANVIS/cand_images/'+ run +'/cand_'+format(ID, '05')+'_'+ field +'_'+ run +'/'
    for cutouts in os.listdir(path_cutout):
        if cutouts.endswith('.png'):
            files.append(path_cutout + cutouts)
    writer = imageio.get_writer(str(path_cutout)  + '_VIDEO.gif', fps =3)
    video_loc = str(path_cutout)  + '_VIDEO.gif'
    for i in files:
        writer.append_data(imageio.imread(i))
    writer.close()

    print('\nDone! Look for your outputs here: /fred/oz100/CANVIS/cand_images/%s' % run)

    print('\n###########################')
    print('#  CANVIS HAS Finished    #')
    print('#  Enjoy and discover!    #')
    print('###########################\n')

    return video_loc


if __name__ == "__main__":

    # Read in arguments
    arguments       = docopt.docopt(__doc__)

    # Mandatory arguments
    field = arguments['<field>']
    ID    = arguments['<ID>']
    ID    = int(ID)
    seed  = arguments['<seed>']
    run   = arguments['<DWF_run>']

    # Optional arguments
    verbose = arguments['--verbose']
    debugmode       = arguments['--debug']
    
    if debugmode:
        print(arguments)

    run_canvis(field,ID,seed,run,verbose=verbose,debugmode=debugmode)
