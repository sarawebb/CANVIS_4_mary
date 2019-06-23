"""
canvis_findingCO.py  -- Input an RA and DEC cut out postage stamps around associated RA and DEC for all available data for given field; Input also the DWF_run, which specifies where CANVIS outputs will be saved. 

Usage: location_canvis [-h] [-v] [--debug] <RA> <DEC> <source_name> <field> <template_image>

Arguments:
        RA (float)
            The RA in degrees.
        DEC (float)
            The DEC in degrees.
        source_name (string)
            The name you would like to call your source in the outputs.
        field (string)
            The DWF field name. 
        template_image (string)
            The template fits file you want to use for the finding chart.

Options:
        -h, --help                              Show this screen
        -v, --verbose                           Show extra information [default: False]     
        --debug                                 Output more for debugging [default: False]

Example:
       python canvis_findingCO.py  -v 246.5000000 -73.0000000 test_cutout ngc6101 REQ16A_20160727_8c29d52-g-20160728T041310_si1_dither_seeing12_join.fits
"""
import docopt
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from astropy.table import Table, Column, join 
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
from astroML.crossmatch import crossmatch_angular 
########################## Manual Input ################################## 


#RA = 246.5000000   
#DEC = -73.0000000
#field = 'ngc6101'
#template = 'REQ16A_20160727_8c29d52-g-20160728T041310_si1_dither_seeing12_join.fits'

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

def create_finding_chart_inputs(RA, DEC,source_name, field, template, verbose=False, debugmode=False):
    print('#############################################')
    print('# CANVIS Finding Chart version has started #')
    print('#############################################')
    print(f'Finding chart image for RA: {RA} DEC: {DEC}')

    temp_path ='/fred/oz100/pipes/DWF_PIPE/TEMPLATES/' + field +'/'+ template
    print(f'Template path being used is: {temp_path}')
    output_path = '/fred/oz100/FINDING_CHARTS/'
    output_files_path = '/fred/oz100/FINDING_CHARTS/'+field+'/'+source_name +'/'
    print(f'Your files will be saved here: {output_files_path}')

    if not os.path.exists(output_files_path):
        os.makedirs(output_files_path, 0o755)
    else:
        pass

    mydic = SortedDict()

    with fits.open(temp_path) as hdu:
        print('## Template opened ##')
        size = 600
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
        hdu.writeto(output_files_path + source_name+'.fits', overwrite = True)
        print('## Fits CUTOUT saved ##')
        plt.axis('off')
        plt.imshow(hdu[0].data, cmap='gray', vmin= -1, vmax =10)
        plt.savefig(output_files_path + source_name+'.png', overwite=True)
        plt.close()
    print('## Starting Source Extracting  ##')

    os.system('sex ' + output_files_path + source_name+'.fits' +' -c default_params.sex -CATALOG_NAME ' + output_files_path + source_name + '_SE.cat -DETECT_THRESH 1.5  -MAG_ZEROPOINT 25 ')
    print('## Source extractor done! ##')
    print('#############################################')
    print('## Cross matching with SkyMapper DR2 ##')
    print('#############################################')
    skymapper_cat_path = '/fred/oz100/SM_cats/'+field+'_SM.csv'
    sm_cat = pd.read_csv(skymapper_cat_path).to_dict(orient='row')
    sm_ra = []
    sm_dec = []
    sm_ngood = []
    sm_class_star = []
    sm_g_psf = []
    sm_g_psf_error =[]
    sm_r_psf = []
    sm_r_psf_error = []
    sm_i_psf = []
    sm_i_psf_error = []
    sm_z_psf = []
    sm_z_psf_error = []

    for i in range(len(sm_cat)):
        sm_ra.append(sm_cat[i]['raj2000'])
        sm_dec.append(sm_cat[i]['dej2000'])
        sm_ngood.append(sm_cat[i]['ngood'])
        sm_class_star.append(sm_cat[i]['class_star'])
        sm_g_psf.append(sm_cat[i]['g_psf'])
        sm_g_psf_error.append(sm_cat[i]['e_g_psf'])
        sm_r_psf.append(sm_cat[i]['r_psf'])
        sm_r_psf_error.append(sm_cat[i]['e_r_psf'])
        sm_i_psf.append(sm_cat[i]['i_psf'])
        sm_i_psf_error.append(sm_cat[i]['e_i_psf'])
        sm_z_psf.append(sm_cat[i]['z_psf'])
        sm_z_psf_error.append(sm_cat[i]['e_z_psf'])
        
    SM_table = Table()
    SM_table['RA'] = sm_ra 
    SM_table['DEC'] = sm_dec
    SM_table['ngood'] = sm_ngood
    SM_table['class_star'] = sm_class_star
    SM_table['g_psf'] = sm_g_psf
    SM_table['e_g_psf'] = sm_g_psf_error
    SM_table['r_psf'] = sm_r_psf
    SM_table['e_r_psf'] = sm_r_psf_error
    SM_table['i_psf'] = sm_i_psf
    SM_table['e_i_psf'] = sm_i_psf_error
    SM_table['z_psf'] = sm_z_psf
    SM_table['e_z_psf'] =sm_z_psf_error
    SM_filtered_ra = []
    SM_filtered_dec = []
    SM_filtered_ngood = []
    SM_filtered_class_star = []
    SM_filtered_g_psf = []
    SM_filtered_g_psf_err = []
    SM_filtered_r_psf = []



    SM_filtered_r_psf_err = []
    SM_filtered_i_psf = []
    SM_filtered_i_psf_err = []
    SM_filtered_z_psf = []
    SM_filtered_z_psf_err = []
                                                                                                                                                        
    for row in SM_table:
        if 0.95 <= row['class_star'] <= 1.0:
            SM_filtered_ra.append(row['RA'])
            SM_filtered_dec.append(row['DEC'])
            SM_filtered_ngood.append(row['ngood'])
            SM_filtered_class_star.append(row['class_star'])
            SM_filtered_g_psf.append(row['g_psf'])
            SM_filtered_g_psf_err.append(row['e_g_psf'])
            SM_filtered_r_psf.append(row['r_psf'])
            SM_filtered_r_psf_err.append(row['e_r_psf'])
            SM_filtered_i_psf.append(row['i_psf'])
            SM_filtered_i_psf_err.append(row['e_i_psf'])
            SM_filtered_z_psf.append(row['z_psf'])
            SM_filtered_z_psf_err.append(row['e_z_psf'])

    
    SM_X = np.empty((len(SM_filtered_ra), 6), dtype=np.float64)
    SM_X[:, 0] = SM_filtered_ra
    SM_X[:, 1] = SM_filtered_dec
    SM_X[:, 2] = SM_filtered_ngood
    SM_X[:, 3] = SM_filtered_g_psf
    SM_X[:, 4] = SM_filtered_r_psf
    SM_X[:, 5] = SM_filtered_i_psf
    
    MAG_APER, MAGERR_APER, MAG_AUTO, MAGERR_AUTO, XPEAK_IMAGE, YPEAK_IMAGE, X_IMAGE, Y_IMAGE, ALPHA_J2000, DELTA_J2000 = np.loadtxt(output_files_path + source_name + '_SE.cat', unpack = True)
    DWF_X = np.empty((len(MAG_APER), 2), dtype=np.float64)
    DWF_X[:, 0] = ALPHA_J2000
    DWF_X[:, 1] = DELTA_J2000

    max_radius = 1/3600 #1 arc second
    dist_between, ind_row = crossmatch_angular(DWF_X, SM_X, max_radius)
    match = ~np.isinf(dist_between)
    
    if len(match) != 0: 
        match_table = Table()
        match_table['matched_true_false'] = match
        match_table['matched_ID'] = ind_row
     
        SM_match_true = []
        SM_row_matched = []
        DWF_g_mags_matched = []
        DWF_g_mags_error_matched = []
        DWF_obs_ra_matched = []
        DWF_obs_dec_matched = []
        for row in match_table:
            if row['matched_true_false'] == True:
                SM_match_true.append(row['matched_true_false'])
                SM_row_matched.append(row['matched_ID'])
        SM_RA = []
        SM_DEC = []
        SM_g_mag = []
        SM_r_mag = []
        SM_i_mag = []
        SM_z_mag = []
        for j in SM_row_matched:

            RA = SM_X[j, 0]
            DEC = SM_X[j, 1]
            g_mag = SM_X[j, 2]
            r_mag = SM_X[j, 3]
            i_mag = SM_X[j, 4]
            z_mag = SM_X[j, 5]
            SM_RA.append(RA)
            SM_DEC.append(DEC)

            SM_g_mag.append(g_mag)
            SM_r_mag.append(r_mag)
            SM_i_mag.append(i_mag)
            SM_z_mag.append(z_mag)

        SM_final_table = Table()
        SM_final_table['RA'] = SM_RA 
        SM_final_table['DEC'] = SM_DEC
        SM_final_table['g_mag'] = SM_g_mag
        SM_final_table['r_mag'] = SM_r_mag
        SM_final_table['i_mag'] = SM_i_mag
        SM_final_table['z_mag'] = SM_z_mag
        
        output= output_files_path + source_name + '_skymapper_star_cat.ascii'
        SM_final_table.write(output, format='ascii', overwrite=True)
        print('#############################################')
        print('# YOUR FINDING CHART INPUTS ARE DONE#')
        print(f'# FIND THEM HERE: {output_files_path}')
        print('#############################################')
if __name__ == "__main__":
    ## read in arguments 
    arguments = docopt.docopt(__doc__, options_first=True)

    # Mandatory arguments
    RA = arguments['<RA>']
    DEC = arguments['<DEC>']
    source_name = arguments['<source_name>']
    field = arguments['<field>']
    print(arguments)
    template_image = arguments['<template_image>']
    
    # Optional arguments
    verbose = arguments['--verbose']
    debugmode   = arguments['--debug']

    if debugmode:
            print(arguments)

    create_finding_chart_inputs(RA,DEC,source_name,field,template_image,verbose=verbose,debugmode=debugmode)
    
