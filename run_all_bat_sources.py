import numpy as np
import pandas as pd
import pickle
import sys

morph_dir = '/Users/ttshimiz/Github/bat-agn-fir-morphology/'
sys.path.append(morph_dir)

import bat_morph_tools as bmt

bat_info = pd.read_csv('/Users/ttshimiz/Github/bat-data/bat_info.csv', index_col=0)
bat_info = bat_info.drop(np.array(['Mrk3', 'NGC3079']))
pacs_info = pd.read_csv(morph_dir+'pacs70_image_info.csv', index_col=0)

names = bat_info.index.values

out_dir = morph_dir+'results/sersic/PACS70/'
models = ['sersic', 'sky']
instrument = 'PACS'
band = 'blue'

result_dict = {}

for n in names:
    
    pra = np.int(pacs_info.loc[n, 'PSF Rotation Angle'])
    result = bmt.run_galfit_single(n, instrument, band, models, region_size=5.0,
                                   psf_rotation=pra, out_dir=out_dir, out_suffix='_sersic')
                                   
    result_dict[n] = result

f = open(morph_dir+'sersic_fit_results_02-29-2016.pkl', 'wb')
pickle.dump(result_dict, f)
f.close()