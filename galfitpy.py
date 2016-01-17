# Functions and script to fit the alphaTau PSF made with scanamorphos with a
# Sersic to the BAT AGN

# Standard scientific imports
import numpy as np
import pandas as pd
import os

# Astropy imports
import astropy.units as u
import astropy.coordinates as coord
import astropy.io.fits as pyf
from astropy import wcs
from astropy.stats import sigma_clipped_stats, gaussian_sigma_to_fwhm, gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel

# Photutils imports
from photutils import detect_sources, segment_properties

# FWHM Sizes for PACS PSF
psf_fwhm = {'blue': 5.299,
            'red': 10.632}

# Pixel Scale
pix_scale = {'blue': 1.4,
             'red': 2.85}


def create_galfit_input_param(filename, src_file, out_file, err_file, psf_file,
                              img_region, conv_box, psf_fine_samp=1.0, bad_pix_mask_file=None,
                              constraint_file=None, zero_pnt=0.0, plate_scale=1.0, display_type='regular',
                              run_option=0):
                              
    fn = open(filename, 'w')
    
    fn.write('================================================================================\n'
              '# IMAGE and GALFIT CONTROL PARAMETERS\n')
              
    fn.write('A) '+src_file+' # Image Map\n')
    fn.write('B) '+out_file+' # Output result file\n')
    fn.write('C) '+err_file+' # Sigma Map\n')
    fn.write('D) '+psf_file+' # PSF Map\n')
    fn.write('E) '+str(psf_fine_samp)+' # PSF fine sampling factor relative to data\n')
    
    if bad_pix_mask_file is not None:
        fn.write('F) '+bad_pix_mask_file+' # Bad pixel mask (FITS image or ASCII coord list)\n')
    else:
        fn.write('F) none # Bad pixel mask (FITS image or ASCII coord list)\n')
    
    if constraint_file is not None:
        fn.write('G) '+constraint_file+' # File with parameter constraints (ASCII file)\n')
    else:
        fn.write('G) none # File with parameter constraints (ASCII file)\n')
    
    fn.write('H) {0:3d} {1:3d} {2:3d} {3:3d} # Image region to fit (xmin xmax ymin ymax)\n'.format(img_region[0], img_region[1], img_region[2], img_region[3]))
    fn.write('I) {0:3d} {1:3d} # Size of the convolution box (x y)\n'.format(conv_box[0], conv_box[1]))
    fn.write('J) '+str(zero_pnt)+' # Magnitude photometric zeropoint\n')
    fn.write('K) {0:0.1f} {1:0.1f} # Plate scale (dx dy)   [arcsec per pixel]\n'.format(plate_scale, plate_scale))
    fn.write('O) '+display_type+' # Display type (regular, curses, both)\n')
    fn.write('P) '+str(run_option)+' # Options: 0=normal run; 1,2=make model/imgblock & quit; 3=seperate model components\n')
    fn.write('\n# Model Components\n')
    
    fn.close() 


def add_sersic_component(input_file, center, int_mag, r_e, sersic_n, axis_ratio, pa,
                         center_fixed=False, int_mag_fixed=False, r_e_fixed=False,
                         sersic_n_fixed=False, axis_ratio_fixed=False, pa_fixed=False,
                         output_option=0):

    fn = open(input_file, 'a')
    
    fn.write('\n# Sersic\n')
    fn.write('0) sersic # Component Type\n')
   
    if center_fixed:
        fn.write('1) {0:0.2f} {1:0.2f} 0 0 # Component center\n'.format(center[0], center[1]))
    else:
        fn.write('1) {0:0.2f} {1:0.2f} 1 1 # Component center\n'.format(center[0], center[1]))

    if int_mag_fixed:
        fn.write('3) {0:0.2f} 0 # Integrated Magnitude\n'.format(int_mag))
    else:
        fn.write('3) {0:0.2f} 1 # Integrated Magnitude\n'.format(int_mag))
    
    if r_e_fixed:
        fn.write('4) {0:0.2f} 0 # Half-light radius\n'.format(r_e))
    else:
        fn.write('4) {0:0.2f} 1 # Half-light radius\n'.format(r_e))
    
    if sersic_n_fixed:
        fn.write('5) {0:0.2f} 0 # Sersic Index n\n'.format(sersic_n))
    else:
        fn.write('5) {0:0.2f} 1 # Sersic Index n\n'.format(sersic_n))
    
    if axis_ratio_fixed:
        fn.write('9) {0:0.2f} 0 # Axis ratio (b/a)\n'.format(axis_ratio))
    else:
        fn.write('9) {0:0.2f} 1 # Axis ratio (b/a)\n'.format(axis_ratio))
    
    if pa_fixed:
        fn.write('10) {0:0.2f} 0 # Position angle\n'.format(pa))
    else:
        fn.write('10) {0:0.2f} 1 # Position angle\n'.format(pa))
    
    fn.write('Z) '+str(output_option)+" #  Output option (0 = resid., 1 = Don't subtract) \n")
    
    fn.close()


def add_gaussian_component(input_file, center, int_mag, fwhm, axis_ratio, pa,
                      center_fixed=False, int_mag_fixed=False, fwhm_fixed=False,
                      axis_ratio_fixed=False, pa_fixed=False, output_option=0):

    fn = open(input_file, 'a')
    
    fn.write('\n# Gaussian\n')
    fn.write('0) gaussian # Component Type\n')
   
    if center_fixed:
        fn.write('1) {0:0.2f} {1:0.2f} 0 0 # Component center\n'.format(center[0], center[1]))
    else:
        fn.write('1) {0:0.2f} {1:0.2f} 1 1 # Component center\n'.format(center[0], center[1]))

    if int_mag_fixed:
        fn.write('3) {0:0.2f} 0 # Integrated Magnitude\n'.format(int_mag))
    else:
        fn.write('3) {0:0.2f} 1 # Integrated Magnitude\n'.format(int_mag))
    
    if fwhm_fixed:
        fn.write('4) {0:0.2f} 0 # FWHM\n'.format(fwhm))
    else:
        fn.write('4) {0:0.2f} 1 # FWHM\n'.format(fwhm))
    
    if axis_ratio_fixed:
        fn.write('9) {0:0.2f} 0 # Axis ratio (b/a)\n'.format(axis_ratio))
    else:
        fn.write('9) {0:0.2f} 1 # Axis ratio (b/a)\n'.format(axis_ratio))
    
    if pa_fixed:
        fn.write('10) {0:0.2f} 0 # Position angle\n'.format(pa))
    else:
        fn.write('10) {0:0.2f} 1 # Position angle\n'.format(pa))
    
    fn.write('Z) '+str(output_option)+" #  Output option (0 = resid., 1 = Don't subtract) \n")
    
    fn.close()

    
def add_psf_component(input_file, center, int_mag, center_fixed=False,
                      int_mag_fixed=False, output_option=0):

    fn = open(input_file, 'a')
    
    fn.write('\n# PSF\n')
    fn.write('0) psf # Component Type\n')
   
    if center_fixed:
        fn.write('1) {0:0.2f} {1:0.2f} 0 0 # Component center\n'.format(center[0], center[1]))
    else:
        fn.write('1) {0:0.2f} {1:0.2f} 1 1 # Component center\n'.format(center[0], center[1]))

    if int_mag_fixed:
        fn.write('3) {0:0.2f} 0 # Integrated Magnitude\n'.format(int_mag))
    else:
        fn.write('3) {0:0.2f} 1 # Integrated Magnitude\n'.format(int_mag))
    
    fn.write('Z) '+str(output_option)+" #  Output option (0 = resid., 1 = Don't subtract) \n")
    
    fn.close()    
    

def add_sky_component(input_file, sky_center, dx=0.0, dy=0.0, sky_center_fixed=False,
                      dx_fixed=True, dy_fixed=True, output_option=0):

    fn = open(input_file, 'a')
    
    fn.write('\n# Sky\n')
    fn.write('0) sky # Component Type\n')
   
    if sky_center_fixed:
        fn.write('1) {0:0.2f} 0 # Sky background at image center\n'.format(sky_center))
    else:
        fn.write('1) {0:0.2f} 1 # Sky background at image center\n'.format(sky_center))
        
    if dx_fixed:
        fn.write('2) {0:0.2f} 0 # sky gradient in x direction\n'.format(dx))
    else:
        fn.write('2) {0:0.2f} 1 # sky gradient in x direction\n'.format(dx))
    
    if dy_fixed:
        fn.write('3) {0:0.2f} 0 # sky gradient in y direction\n'.format(dy))
    else:
        fn.write('3) {0:0.2f} 1 # sky gradient in y direction\n'.format(dy))
    
    fn.write('Z) '+str(output_option)+" #  Output option (0 = resid., 1 = Don't subtract) \n")
    
    fn.close()    
    
 
# Function to estimate the global sky background
def estimate_sky(image, clip=3.0):

    # Create the NaN mask
    nanmask = np.isnan(image)
    
    # Estimate of the global sky background from sigma-clipping
    im_mean, im_med, im_std = sigma_clipped_stats(image, sigma=clip, mask=nanmask)
    
    return im_med, im_std


# Function to find the source from segmentation map
def find_source(img, thresh, im_wcs, coord_src, rdist):

    segm_img = detect_sources(img, thresh, npixels=5)
    props = segment_properties(img, segm_img, wcs=im_wcs)
    ra_src = [p.ra_icrs_centroid for p in props]
    dec_src = [p.dec_icrs_centroid for p in props]

    coord_src_img = coord.SkyCoord(ra_src, dec_src, frame='icrs')
    dist = coord_src.separation(coord_src_img).to(u.arcsec)
    ind_closest = np.argmin(dist)
    cdist = np.min(dist).value

    if cdist < rdist:
        return props[ind_closest]
    else:
        return None   


# Function that converts the source properties from photutils to initial parameters
# for GALFIT
def get_galfit_param(img, src, waveband, region_size=10., model='sersic'):

    center = [src.xcentroid.value, src.ycentroid.value]
    if model != 'psf':
        axis_ratio = 1./src.elongation.value
        pa = src.orientation.to(u.deg).value + 90
        fwhm_major = src.semimajor_axis_sigma.value*2.355
        fwhm_minor = src.semiminor_axis_sigma.value*2.355
        fwhm_circle = np.sqrt((fwhm_major**2 + fwhm_minor**2)/2.)
    
        if fwhm_circle > psf_fwhm[waveband]/pix_scale[waveband]:
            fwhm_int = np.sqrt(fwhm_circle**2 - (psf_fwhm[waveband]/pix_scale[waveband])**2)
        else:
            fwhm_int = 1.0
    
        if model=='gaussian':
            src_size = fwhm_int
        elif model=='sersic':
            src_size = fwhm_int/2.0
        
        if src_size < 1.0:
            src_size = 1.0
        
        box_size = np.max([region_size*fwhm_int, 75.0])
    else:
        box_size = 75.0

    region = np.array([center[0]-box_size/2., center[0]+box_size/2.,
                       center[1]-box_size/2., center[1]+box_size/2.],
                      dtype=np.int)
    if region[0] < 1:
		region[0] = 1
    if region[2] < 1:
        region[2] = 1
    if region[1] > img.shape[1]:
        region[1] = img.shape[1]
    if region[3] > img.shape[0]:
        region[3] = img.shape[0]
	
    if model == 'psf':
	    return center, region
    else:
        return center, region, src_size, axis_ratio, pa

# Function to create a bad pixel mask. Outputs a file that is a list
# of pixel coordinates that are either NaN in the signal map or 
# 0/infinite in the error map.
 def create_bad_pixel_mask(img, err_img, filename='bad_pixels.txt'):

    nanmask = np.isnan(img)
    errmask = (err_img <= 0) | (~np.isfinite(err_img))
    fn = open(filename, 'w')
    
    for x in range(img.shape[1]):
        for y in range(img.shape[0]):
            if (nanmask[y, x] | errmask[y, x]):
                fn.write('{0:d} {1:d}\n'.format(x+1, y+1))
    
    fn.close()
    
    
def run_galfit(input_file):

    os.system('~/src/galfit3.0/galfit '+input_file)