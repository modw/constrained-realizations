import healpy as hp
import numpy as np
import mm_constrained_realizations as mmcr

import tempfile
import sys
import os

# params
nside = 64
fwhm = (160/60)*np.pi/180 # from ISW fits file

# maps

weights_map = hp.read_map('input_files/full_weights_map_08_0064.fits',
                          verbose=False)
cl = np.loadtxt('input_files/Cl_ISW_scalar_PlanckPR2ISWpaper.dat')
isw, stdev, mask = hp.read_map('input_files/COM_CompMap_ISW_0064_R2.00.fits',
                               hdu=1,field={0,1,2},verbose=False)

# converting units to muK
tcmb = 2.7255 # from Planck params paper
isw[mask == 1] *= 2.7255*10**6
stdev[mask==1] *= 2.7255*10**6
# noise_cov
noise_cov = stdev**2
noise_cov[mask==0] = hp.UNSEEN

### initializing CR
cr = mmcr.ConstrainedRealizations(nside,weights_map,mask=mask)
# setting covariances
cr.set_signal_cov(cl,fwhm)
cr.set_noise_cov(noise_cov)

# cooling schedule params
lamb_0 = 1e4
tp = 1e-8
cr.set_cooling_schedule(lamb_0,tp)

# number of realizations
n = int(sys.argv[1])

# create temp folder
fdir = tempfile.mkdtemp(prefix='ISW_cr_PPR2_hdu1_',
                        dir=os.getcwd())
# save maps in temp_folder
for i in range(n):
    #create file name
    fd, fname = tempfile.mkstemp(prefix='ISW_cr_PPR2_hdu1_', dir=fdir,
                                 suffix='.npy')
    flm = cr.gen_constrained_realization()
    np.save(fname,flm)
