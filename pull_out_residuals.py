#Script to pull out the residual image for all of the GALFIT results files

import astropy.io.fits as pyf
import numpy as np
import pandas as pd
import pickle
import sys
import os

morph_dir = '/Users/ttshimiz/Github/bat-agn-fir-morphology/'
sys.path.append(morph_dir)

bat_info = pd.read_csv('/Users/ttshimiz/Github/bat-data/bat_info.csv', index_col=0)
bat_info = bat_info.drop(np.array(['Mrk3', 'NGC3079']))
pacs_info = pd.read_csv(morph_dir+'pacs70_image_info.csv', index_col=0)

names = bat_info.index.values
#names = ['Mrk79']
out_dir = morph_dir+'results/sersic/PACS70/'

for n in names:
    
    fn = out_dir+n+'_pacs70_sersic.fits'

    if os.path.exists(fn):  
        print 'Pulling out residuals for ', n
        residuals = pyf.getdata(fn, ext = 3)
        header = pyf.getheader(fn, ext = 2)
        header_image = pyf.getheader(fn, ext = 1)
        fit_region = header['FITSECT']
        xshift = int(fit_region.split(',')[0].split(':')[0][1:])
        yshift = int(fit_region.split(',')[1].split(':')[0])
        header['BUNIT'] = header_image['BUNIT']
        header['PFOV'] = header_image['PFOV']
        header['CTYPE1'] = header_image['CTYPE1']
        header['CTYPE2'] = header_image['CTYPE2']
        header['EQUINOX'] = header_image['EQUINOX']
        header['CRPIX1'] = header_image['CRPIX1'] - xshift+1
        header['CRPIX2'] = header_image['CRPIX2'] - yshift+1
        header['CRVAL1'] = header_image['CRVAL1']
        header['CRVAL2'] = header_image['CRVAL2']
        header['CDELT1'] = header_image['CDELT1']
        header['CDELT2'] = header_image['CDELT2']
        header['CROTA2'] = header_image['CROTA2']
        hdu = pyf.PrimaryHDU(data = residuals, header = header)
        fn_split = fn.split('.')
        new_fn = fn_split[0]+'_residuals.fits'
        hdu.writeto(new_fn, clobber = True)
    
