# Functions and script to fit the alphaTau PSF made with scanamorphos with a
# Sersic to the BAT AGN

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

