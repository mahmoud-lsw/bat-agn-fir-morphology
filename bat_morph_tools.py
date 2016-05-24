import numpy as np
import pandas as pd
import astropy.coordinates as coord
import astropy.wcs as wcs
import astropy.io.fits as fits
import astropy.units as u
import os
import sys

sys.path.append('../galfitpy/')
import galfitpy

import matplotlib.pyplot as plt
import aplpy

wavebands = {'PACS': {'blue': '70',
                      'red': '160'},
             'SPIRE': {'blue': '250',
                       'green': '350',
                       'red': '500'}}

pix_scale = {'PACS': {'blue': 1.4,
                      'red': 2.85},
             'SPIRE': {'blue': 4.5,
                       'green': 6.25,
                       'red': 9.0}}

psf_fwhm = {'PACS': {'blue': 6.26,
                     'red': 12.30},
            'SPIRE': {'blue': 17.6,
                      'green': 23.9,
                      'red': 35.2}}
                       

def get_bat_coords(name):

    bat_info = pd.read_csv('/Users/ttshimiz/Github/bat-data/bat_info.csv',
                           index_col=0)
    cs= coord.SkyCoord(bat_info.loc[name, 'RA_deg']*u.deg,
                       bat_info.loc[name, 'DEC_deg']*u.deg)

    return cs


def get_bat_image(name, instrument, band):

    dir = '/Users/ttshimiz/Dropbox/Herschel_Images/'+instrument+'/'
    
    wave = wavebands[instrument][band]
    
    src_file = dir+name+'_scanamorphos_'+instrument+wave+'_signal.fits'
    err_file = dir+name+'_scanamorphos_'+instrument+wave+'_error.fits'
    
    hdu_src = fits.open(src_file)
    hdu_err = fits.open(err_file)
    
    src_img = hdu_src[0].data
    err_img = hdu_err[0].data
    wcs_img = wcs.WCS(hdu_src[0].header)
    
    return src_img, err_img, wcs_img
    
def get_best_fit_params(output_file, models):

    # Setup dictionary to hold parameters
    params = {}
    if os.path.exists(output_file):
        # Best fit parameters are always in the 3rd frame's header
        hdu = fits.open(output_file)
        header = hdu[2].header       
    
        for i,mod in enumerate(models):
            params[mod] = {}
        
            # Sersic Component
            if mod == 'sersic':
                ncomp = str(i+1)
                xstr = header[ncomp+'_XC']        # X-pixel of center
                ystr = header[ncomp+'_YC']        # Y-pixel of center 
                magstr = header[ncomp+'_MAG']     # Integrated magnitude
                restr = header[ncomp+'_RE']       # Half-light radius
                nstr = header[ncomp+'_N']         # Sersic index
                arstr = header[ncomp+'_AR']       # Axis ratio
                pastr = header[ncomp+'_PA']       # Position angle
            
                # Extract the best fit parameter value and its error
                # If the parameter was fixed then the error will be 0
                params[mod]['center_x'] = parse_param_str(xstr)
                params[mod]['center_y'] = parse_param_str(ystr)
                params[mod]['mag'] = parse_param_str(magstr)
                params[mod]['re'] = parse_param_str(restr)
                params[mod]['n'] = parse_param_str(nstr)
                params[mod]['axis_ratio'] = parse_param_str(arstr)
                params[mod]['pa'] = parse_param_str(pastr)
        
            # PSF component    
            elif mod == 'psf':
                ncomp = str(i+1)
                xstr = header[ncomp+'_XC']
                ystr = header[ncomp+'_YC']
                magstr = header[ncomp+'_MAG']
   
                params[mod]['center_x'] = parse_param_str(xstr)
                params[mod]['center_y'] = parse_param_str(ystr)
                params[mod]['mag'] = parse_param_str(magstr)
        
            # Gaussian component    
            elif mod == 'gauss':
        
                ncomp = str(i+1)
                xstr = header[ncomp+'_XC']
                ystr = header[ncomp+'_YC']
                magstr = header[ncomp+'_MAG']
                fwhmstr = header[ncomp+'_FWHM']
                arstr = header[ncomp+'_AR']
                pastr = header[ncomp+'_PA']
            
                params[mod]['center_x'] = parse_param_str(xstr)
                params[mod]['center_y'] = parse_param_str(ystr)
                params[mod]['mag'] = parse_param_str(magstr)
                params[mod]['fwhm'] = parse_param_str(restr)
                params[mod]['axis_ratio'] = parse_param_str(arstr)
                params[mod]['pa'] = parse_param_str(pastr)
        
            # Sky component
            elif mod == 'sky':
            
                ncomp = str(i+1)
                sky_center_str = header[ncomp+'_SKY']
                skydxstr = header[ncomp+'_DSDX']
                skydystr = header[ncomp+'_DSDY']
            
                params[mod]['sky_center'] = parse_param_str(sky_center_str)
                params[mod]['sky_dx'] = parse_param_str(skydxstr)
                params[mod]['sky_dy'] = parse_param_str(skydystr)
            
    
        # General fitting parameters
        params['chi2'] = header['CHISQ']      # Total chi-square
        params['ndof'] = header['NDOF']       # Number of degrees of freedom
        params['nfree'] = header['NFREE']     # Number of free parameters
        params['nfix'] = header['NFIX']       # Number of fixed parameters
        params['chi2nu'] = header['CHI2NU']   # Reduced chi-square (chi2/ndof)
    
    return params
    
    

def parse_param_str(param_str):
    
    # Function to extract the best fit parameter value and its error.
    # A fixed parameter is enclosed in square brackets and will be assigned
    # a value of 0.
    value_str = param_str.split()[0]
    if value_str[0] == '[':
    	value = np.float(value_str[1:-1])
        err = 0.0
    elif value_str[0] == '*':
        value = np.float(value_str[1:-1])
        err = np.inf
    else:
        value = np.float(value_str)
        err = np.float(param_str.split()[2])
    
    return value, err
    
def estimate_mag(name, instr, wave):

    band = {'PACS': {'blue': 'PACS70',
                     'red': 'PACS160'},
            'SPIRE': {'blue': 'PSW',
                      'green': 'PMW',
                      'red': 'PLW'}}

    bat_herschel = pd.read_csv('/Users/ttshimiz/Github/bat-data/bat_herschel.csv',
                               index_col=0)
    flux = bat_herschel.loc[name, band[instr][wave]]
    if flux == 0:
       flux = bat_herschel.loc[name, band[instr][wave]+'_err'] 
    mag = -2.5*np.log10(flux)
   
    return mag
   
         
def run_galfit_single(name, instr, wave, models, conv_box=[200,200], region_size=10.,
                      thresh_level=3, search_radius=2, psf_rotation=180,
                      out_suffix='_out', out_dir=''):
    
    fit_notes = []
    
    # Get the image and error data
    src_img, err_img, wcs_img = get_bat_image(name, instr, wave)

    # Get the coordinates of the source
    cs = get_bat_coords(name)
    
    # Estimate the background
    im_med, im_std = galfitpy.estimate_sky(src_img)
    thresh = im_med + thresh_level*im_std
    
    # Find and estimate parameters of the source
    rdist = search_radius*psf_fwhm[instr][wave]
    source = galfitpy.find_source(src_img, thresh, wcs_img, cs, rdist)
    
    # Convert to parameters that GALFIT needs
    if source is not None:
        center, region, src_size, axis_ratio, pa = galfitpy.get_galfit_param(src_img, source, wave,
                                                                             region_size=region_size)                                                                            
    else:
        pix_src = cs.to_pixel(wcs_img, origin=0)
        center = np.array([pix_src[0], pix_src[1]])
        src_size = 1.0
        axis_ratio = 1.0
        pa = 0.0
        
        region = np.array([center[0]-75./2., center[0]+75./2.,
                           center[1]-75./2., center[1]+75./2.],
                           dtype=np.int)
        fit_notes.append('Source not detected.') 

    if (np.any(np.array(models)=='sersic') & (src_size > 1.0)):
        src_size = src_size/2.0
    
    # Initialize the GALFIT input file
    dir = '/Users/ttshimiz/Dropbox/Herschel_Images/'+instr+'/'
    
    wavelength = wavebands[instr][wave]
    
    # Create bad pixel file
    bp_file = out_dir+name+'_'+instr.lower()+wavelength+'_bad_pixels.txt'
    galfitpy.create_bad_pixel_mask(src_img, err_img, filename=bp_file)
    
    # Filenames of the image, error image, psf, and output file
    src_file = dir+name+'_scanamorphos_'+instr.lower()+wavelength+'_signal.fits'
    err_file = dir+name+'_scanamorphos_'+instr.lower()+wavelength+'_error.fits'
    psf_file = 'alphaTaupsf/alphaTau20'+wave+'+'+str(np.int(psf_rotation))+'_scanamorphos.fits'
    out_file = out_dir+name+'_'+instr.lower()+wavelength+out_suffix+'.fits'
    input_file = out_dir+name+'_'+instr.lower()+wavelength+out_suffix+'.input'
    
    plate_scale = pix_scale[instr][wave]
    galfitpy.create_galfit_input_param(input_file, src_file, out_file,
                                       err_file, psf_file, region, conv_box,
                                       bad_pix_mask_file=bp_file, plate_scale=plate_scale) 
    
    # Get an estimate of the magnitude from aperture photometry measurements
    mag = estimate_mag(name, instr, wave)
                                    
    # Add in the model components
    if (models == 'psf'):
        galfitpy.add_psf_component(input_file, center, mag)
    elif (models == 'sersic'):
        galfirpy.add_sersic_component(input_file, center, mag, src_size, 1.0, axis_ratio, pa)
    elif (models == 'gauss'):
        galfitpy.add_sersic_component(input_file, center, mag, src_size, axis_ratio, pa)
    else:
    
        for i in models:
            if (i == 'psf'):
                galfitpy.add_psf_component(input_file, center, mag+2.5)
            elif (i == 'sersic'):
                galfitpy.add_sersic_component(input_file, center, mag, src_size, 1.0, axis_ratio, pa)
            elif (i == 'gauss'):
                galfitpy.add_gaussian_component(input_file, center, mag, src_size, axis_ratio, pa)
            elif (i == 'sky'):
                galfitpy.add_sky_component(input_file, im_med)    

    # Run GALFIT
    galfitpy.run_galfit(input_file)
    
    #Retrieve the best-fit parameters
    result = get_best_fit_params(out_file, models)
    if len(result) > 0:
		if ((models == 'psf') | (models == 'sersic') | (models == 'gauss')):
			keys = result[models].keys()
			test_fit = True
			for k in keys:
				if np.isinf(result[models][k][1]):
					test_fit = False
					break
		else:
			test_fit = True
			for m in models:
				keys = result[m].keys()
				for k in keys:
					if np.isinf(result[m][k][1]):
						test_fit = False
						break
				if not test_fit:
					break
	
		if not test_fit:
			fit_notes.append('At least one problematic parameter.')
    else:
        fit_notes.append('Fit failed')
        
    result['fit_notes'] = fit_notes            
            
    return result


def plot_galfit_results(out_file, name):

    fig = plt.figure(figsize=(24, 6))
    
    data_fig = aplpy.FITSFigure(out_file, hdu=1, subplot=(1,3,1), figure=fig)
    model_fig = aplpy.FITSFigure(out_file, hdu=2, subplot=(1,3,2), figure=fig)
    resid_fig = aplpy.FITSFigure(out_file, hdu=3, subplot=(1,3,3), figure=fig)
    
    data_fig.set_theme('publication')
    model_fig.set_theme('publication')
    resid_fig.set_theme('publication')
    
    data_fig.show_colorscale(cmap='viridis', stretch='log', vmin=0, vmid=-1e-3)
    model_fig.show_colorscale(cmap='viridis', stretch='log', vmin=0, vmid=-1e-3)
    resid_fig.show_colorscale(cmap='viridis', stretch='linear')
    
    data_fig.show_colorbar()
    model_fig.show_colorbar()
    resid_fig.show_colorbar()
    
    data_fig.add_label(0.05, 0.95, name, color='white', relative=True, size=20,
                       horizontalalignment='left', verticalalignment='top')
    
    return fig


def calc_sersic_radius_perc(n, frac):
    """Calculates the normalized radius of a Sersic profile that contains
       fractional flux, 'frac'
    """
    
    # First find kappa such that Gamma(2n) = 2*gamma(2n, kappa)
    # Note that scipy.special.gammainc(a, x) is normalized by 1/Gamma(a)
    g2n = spec.gamma(2*n)
    f1 = lambda k: 2*spec.gammainc(2*n, k)*g2n - g2n
    
    # Use MacArthur+2003 for guess of kappa
    if n > 0.36:
        kappa0 = (2*n - 1./3. + 4./(405*n) + 46./(25515*n**2) + 131./(1148175*n**3) -
                2194697./(30690717750.*n**4))
    else:
        a0 = 0.01945
        a1 = -0.8902
        a2 = 10.95
        a3 = -19.67
        a4 = 13.43
        kappa0 = a0 + a1*n + a2*n**2 + a3*n**3 + a4*n**4
    
    kappa = opt.newton(f1, kappa0, maxiter=1000)
    
    # Solve for the normalized radius = r/re
    # Use Brent's method for root-finding
    # Need to bracket the real value by finding where the sign changes
    f2 = lambda x: spec.gammainc(2*n, kappa*x**(1/n)) - frac
    rtest = np.arange(0, 100., 0.01)
    ftest = f2(rtest)
    a = rtest[ftest < 0][-1] 
    b = rtest[ftest > 0][0]
    rp = opt.brentq(f2, a, b, maxiter=1000)
    
    return rp
    
                                                    