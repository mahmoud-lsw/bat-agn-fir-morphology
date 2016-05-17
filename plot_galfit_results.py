import numpy as np
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt

morph_dir = '/Users/ttshimiz/Github/bat-agn-fir-morphology/'
sys.path.append(morph_dir)

import bat_morph_tools as bmt

bat_info = pd.read_csv('/Users/ttshimiz/Github/bat-data/bat_info.csv', index_col=0)
bat_info = bat_info.drop(np.array(['Mrk3', 'NGC3079']))
pacs_info = pd.read_csv(morph_dir+'pacs70_image_info.csv', index_col=0)

names = bat_info.index.values
#names = ['PG2304+042', 'NGC7582', 'NGC7469', 'NGC6860',
#         'NGC6240', 'NGC5728', 'NGC5506', 'NGC4151', 'NGC4051',
#         'NGC3783', 'NGC3281', 'NGC3081', 'NGC2992', 'NGC1365', 'Mrk841',
#         'Mrk348', 'Mrk335', 'Mrk290', 'IISZ010', 'IC5063', 'CenA']
out_dir = morph_dir+'results/sersic/PACS70/'

for n in names:
    
    out_file = out_dir+n+'_pacs70_sersic.fits'
    pra = np.int(pacs_info.loc[n, 'PSF Rotation Angle'])
    if (os.path.exists(out_file)) & (pra != 180):
        fig = bmt.plot_galfit_results(out_file, n)
        fig.savefig(out_dir+n+'_pacs70_sersic.pdf', bbox_inches='tight')
        plt.close(fig)