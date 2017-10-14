import healpy as hp
import numpy as np
import mm_constrained_realizations as mmcr

import tempfile
import sys
import os

# params
nside = 64

# maps

weights_map = hp.read_map('Commander-Ruler/input_files/full_weights_map_08_0064.fits', verbose=False)
cl = np.load('Commander-Ruler/input_files/cls_PlanckPR2_TT_lowp_lensing_lensed.npy')
tmap = hp.read_map('Commander-Ruler/input_files/commander_t_map.fits', verbose=False)
noise_cov = hp.read_map('Commander-Ruler/input_files/commander_noise_cov.fits', verbose=False)
mask = hp.read_map('Commander-Ruler/input_files/commander_mask.fits', verbose=False)

### initializing CR
cr = mmcr.ConstrainedRealizations(nside,weights_map,mask=mask)
# setting covariances
cr.set_signal_cov(cl)
cr.set_noise_cov(noise_cov)

# cooling schedule params
lamb_0 = 1e10
tp = 1e-6
cr.set_cooling_schedule(lamb_0,tp)

# number of realizations
n = int(sys.argv[1])

# create temp folder
fdir = tempfile.mkdtemp(prefix='commander_ruler_nside64_no_smoothing_CRs_folder_',
                        dir=os.getcwd())
# save maps in temp_folder
for i in range(n):
    #create file name
    fd, fname = tempfile.mkstemp(prefix='commander_ruler_nside64_no_smoothing_CRs_map_', dir=fdir,
                                 suffix='.npy')
    flm = cr.gen_constrained_realization()
    np.save(fname,flm)
