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

names = bat_info.index.values

out_dir = morph_dir+'results/sersic/PACS70/'

for n in names:
    
    out_file = out_dir+n+'_pacs70_sersic.fits'
    if os.path.exists(out_file):
        fig = bmt.plot_galfit_results(out_file, n)
        fig.savefig(out_dir+n+'pacs70_sersic.pdf', bbox_inches='tight')
        plt.close(fig)