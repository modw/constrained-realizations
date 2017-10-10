import numpy as np
import healpy as hp


class Parameters:
    def __init__(self, nside, lmax_factor):
        """Set base parameters"""
        self.nside = nside
        self.npix = hp.nside2npix(nside)
        self.lmax = np.int0(lmax_factor*nside)
        self.pix_area = hp.pixelfunc.nside2pixarea(nside)


class Mask:
    """Instantiate mask, get good_pix and bad_pix and degrade with 0.9 criteria
    if nside_out is specified."""
    def __init__(self, mask, nside_out=None):
        msk = np.copy(mask)
        if nside_out is not None:
            msk = hp.ud_grade(msk, nside_out)
            msk[msk >= 0.9] = 1
            msk[msk < 0.9] = 0
        self.good_pix = msk == 1
        self.bad_pix = msk == 0
        self.mask = msk


class NoiseCov:
    """Noise covariance from map and mask, assigns tau as an attribute"""
    def __init__(self, noise_cov, mask):
        self.tau = np.min(noise_cov[mask == 1])
        self.t_cov = self.tau*np.ones(hp.get_map_size(noise_cov))
        self.noise_cov = noise_cov


class SignalCov:
    """ Set signal Cls and calculate Sl from Cls and lmax """
    def __init__(self, cl, lmax):
        self.cl = cl[:lmax+1]
        self.cl_inv = np.zeros(lmax+1)
        self.cl_inv[2:] = 1/self.cl[2:]
        l_list, m_list = hp.sphtfunc.Alm.getlm(lmax)
        self.signal_cov = cl[l_list]
        self.signal_cov_inv = np.zeros(len(l_list))
        self.signal_cov_inv[l_list > 1] = 1/self.signal_cov[l_list > 1]


class CoolingSchedule:
    def __init__(self, lamb_0, eta):
        self.set_lambda_schedule(lamb_0, eta)

    def set_lambda_schedule(self, lamb_0, eta):
        """
        Get whole lambda> list by multiplying lamb_0 by powers of eta.
        """
        self.lamb_0 = lamb_0
        self.eta = eta
        self.lamb_list = np.array([])

        # constructing list
        lamb = lamb_0
        while lamb > 1:
            self.lamb_list = np.append(self.lamb_list, lamb)
            lamb *= eta
        self.lamb_list = np.append(self.lamb_list, 1)

    def set_precision_schedule(self, target_precision):
        """
        Function to vary target precision as lambda changes.
        Split eps in three parts by factors of 10.
        """
        self.tp = target_precision
        self.eps_list = 10**(2)*target_precision*np.ones(len(self.lamb_list))
        self.eps_list[len(self.lamb_list)//3:] = 10**1 * target_precision
        self.eps_list[2*len(self.lamb_list)//3:] = target_precision


class Delta:
    """Generate random field delta of which the wiener filter is f.
    f is the fluctuation field to be added to the wiener filter of data
    to obtain one constrained realization."""
    def __init__(self, noise_cov, cl_inv, mask, pix_area):
        self.sqrt_ncov_inv = 1/np.sqrt(noise_cov[mask == 1])
        self.cl_g2s = (pix_area**2)*cl_inv

    def gen_delta(self, good_pix, bad_pix, nside, npix):
        g1n = self.sqrt_ncov_inv*np.random.normal(size=len(self.sqrt_ncov_inv))
        g2s = hp.synfast(self.cl_g2s, nside, len(self.cl_g2s)-1, verbose=False)
        delta = np.zeros(npix)
        delta[good_pix] = g1n + g2s[good_pix]
        delta[bad_pix] = g2s[bad_pix]
        return delta
