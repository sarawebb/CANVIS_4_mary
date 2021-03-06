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

ID = 21847
field = 'test'
seed ='NOAO5'
run ='june18'

parser = ArgumentParser()
parser.add_argument('--ID', required=False, 
                        help='ID of the object you want to use CANVIS on.',
                        dest='ID', nargs='+')
			
parser.add_argument('--Field', required=True, 
                        help='field the object is in.',
                        dest='field', nargs='+')
			
parser.add_argument('--run', required=False, 
                        help='What DWF run?',
                        dest='run', nargs='+')

parser.add_argument('--seed', required=False, 
                        help='What the seed mary run?',
                        dest='seed', nargs='+')


translist_path = '/home/fstars/MARY4OZ/transients/transients_coo_'+field+'_'+seed+'.txt'

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

ra_list = []
dec_list = []
field_list = []
ccd_num_list = []

with open(translist_path) as f:
    for line in f:
        line = line.split()
        if(int(line[0]) in objects):
            ra_list.append(str(line[1]))
            dec_list.append(str(line[2]))
            field_list.append(str(line[3]))            
            ccd_num_list.append(str(line[6]))   

for objid, ra, dec, field, ccd_num in zip(objects, ra_list, dec_list, field_list, ccd_num_list):
    print(ra, dec, field, ccd_num)


    path ='/fred/oz100/pipes/DWF_PIPE/MARY_WORK/'+field+'_*_*_*/ccd' + ccd_num+'/images_resampled/sci_' + ccd_num+'.resamp.fits'
    print(path)
    fitsfileslist = glob.glob(path)
    print(fitsfileslist)
    #print(fitsfileslist)
    mydic = SortedDict()
    for path in fitsfileslist:
    	for fitsim in os.listdir(path):
            with fits.open(path+fitsim) as hdu:
                size = 200 
                w = WCS(hdu[0].header)
                #print(w)
                head = hdu[0].header
                xlim=head['NAXIS1']
                ylim=head['NAXIS2']
                #print(xlim, ylim)
                pixcrd_im = np.array([[xlim, ylim]], np.float_)
                world_im = w.wcs_pix2world(pixcrd_im, 1)
                pixx_im, pixy_im = world_im[0][0], world_im[0][1]
                corners=w.calc_footprint()
                #print(corners)
                corner_1 = corners[0]
                corner_2 = corners[1]
                corner_3 = corners[2]
                corner_4 = corners[3]
                differnce = corner_1 - corner_2 
                #print(differnce)
                #print(corner_1)
                
                    
                pixcrd = np.array([[RA_test, DEC_test]], np.float_)
                #print(pixcrd)
                worldpix = w.wcs_world2pix(pixcrd, 1)
                #print(worldpix)
                pixx, pixy = worldpix[0][0], worldpix[0][1]
                #print(pixx, pixy)
    
                #print(pixx, pixy)
                print('before loop')
        
                #print(RA_test, DEC_test)
                #print(pixx_im,pixy_im)
        
                #print(corner_1[0], corner_4[0])
                #print(corner_1[0],  RA_test,corner_2[0])
                #print(corner_1[1],  DEC_test,corner_3[1])
                if  corner_1[0] <= RA_test <=corner_2[0] and corner_1[1] >= DEC_test >= corner_3[1]:
                    mydic[date] =[path, pixx, pixy]
    for i, (key, (path, pixx, pixy)) in enumerate(mydic.items()):
        path_cand = '/fred/oz100/CANVIS/cand_images/'+ run + '/cand_'+format(objid, '05')+'_'+ field +'_'+ run +'/'
        path_cutout = '/fred/oz100/CANVIS/cand_images/'+ run +'/cand_'+format(objid, '05')+'_'+ field +'_'+ run +'/cand_'+format(objid, '05')+'_'+run+'_cutout_'+format(i, '03')+'.fits'  
        if not os.path.exists(path_cand):
            os.makedirs(path_cand, 0o755)
        else:
            pass
        size = 200
        with fits.open(path) as hdu:
            nom_data = (hdu[0].data - np.min(hdu[0].data))/(np.max(hdu[0].data)-np.min(hdu[0].data))
            print(' FOUND ON THIS CCD ' +fitsim )
            cutout = Cutout2D(hdu[0].data, (pixx, pixy), size, wcs= w)
            hdu[0].data = cutout.data
            hdu[0].header['CRPIX1'] = cutout.wcs.wcs.crpix[0]
            hdu[0].header['CRPIX2'] = cutout.wcs.wcs.crpix[1]
            hdu.writeto(path_cutout+ fitsim+'_cutout.fits', overwrite = True)
            plt.axis('off')
            plt.imshow(hdu[0].data, cmap='gray')
            plt.colorbar()
            plt.savefig(path_cutout+'cutout_'+fitsim+'.png')
            plt.close()
	
	
	
    
    print('Hello')
    files = []	
    for cutouts in os.listdir(path_cutout):
    	if cutouts.endswith('.png'):
            files.append(path_cutout + cutouts)
    writer = imageio.get_writer(path_cutout + i + '_VIDEO.gif', fps =3)
    for i in files:
    	writer.append_data(imageio.imread(i))
    writer.close()
	

