# messenger method (mm) constrained realizations

import numpy as np
import healpy as hp
import support_classes as sc


class ConstrainedRealizations:
    def __init__(self, nside, weights_map, lmax_factor=1.5, mask=None):
        """Constrained Realizations has to be initialized with parameters and
        the desired mask."""
        self.params = sc.Parameters(nside, lmax_factor)
        self.weights_map = weights_map
        if mask is not None:
            self.mask = sc.Mask(mask, nside)
        else:
            self.mask = sc.Mask(np.ones(self.params.npix))

    # set up methods
    def set_noise_cov(self, noise_cov):
        """Set noise covariance matrix and calculate quantities for messenger
        loop."""
        self.check_nside(noise_cov)
        self.noise_cov = sc.NoiseCov(noise_cov, self.mask.mask)
        try:
            self.set_delta()
        except:
            pass

    def set_signal_cov(self, cl, fwhm = None):
        """Set signal cov from Cls and assignt quantities to .signal_cov
        atribute"""
        if fwhm is not None:
            bl = hp.gauss_beam(fwhm,lmax=len(cl)-1)
            self.signal_cov = sc.SignalCov(bl**2 * cl, self.params.lmax)
        else:
            self.signal_cov = sc.SignalCov(cl, self.params.lmax)
        try:
            self.set_delta()
        except:
            pass

    def set_delta(self):
        """Set delta to generate random filed. Use gen_delta method to generate
        delta"""
        self.delta = sc.Delta(self.noise_cov.noise_cov, self.signal_cov.cl_inv,
                              self.mask.mask, self.params.pix_area)

    def set_cooling_schedule(self, lamb_0, target_precision, eta=3/4):
        """Set initial parameters to obtain cooling schedule under atribute .cs.
        """
        self.cs = sc.CoolingSchedule(lamb_0, eta)
        self.cs.set_precision_schedule(target_precision)

    # cooling methods/functions
    def gen_delta(self):
        """Generate one fluctuation field delta making use of already set
        parameters and covariances."""
        delta = self.delta.gen_delta(self.mask.good_pix, self.mask.bad_pix,
                                     self.params.nside, self.params.npix)
        return delta

    def do_transform(self, delta, tlm, lamb):
        """Do one iteration of the basis transform for the wiener filter"""
        sl, tau, noise_cov = \
            self.signal_cov.signal_cov, self.noise_cov.tau, \
            self.noise_cov.noise_cov
        good_pix, bad_pix = self.mask.good_pix, self.mask.bad_pix
        weights_map = self.weights_map
        pix_area, lmax, nside = self.params.pix_area, self.params.lmax, \
            self.params.nside
        flm, tlm = \
            field_trnsfrm(delta, tlm, lamb, sl, tau, noise_cov, good_pix,
                          bad_pix, weights_map, pix_area, lmax, nside)
        return flm, tlm

    def solve_flm(self, tlm, delta, lamb, target_precision):
        """Solve for flm for a specific delta, lambda and target precision"""
        flm = np.ones(len(tlm))
        while True:
            flm_i = np.copy(flm)
            flm, tlm = self.do_transform(delta, tlm, lamb)
            conv = np.linalg.norm(flm-flm_i, ord=2)\
                / np.linalg.norm(flm_i, ord=2)
            if conv < target_precision:
                return flm, tlm

    def gen_constrained_realization(self, delta_fix=None):
        """Generate one fluctation field, which added to a wiener filtered map
        gives a constrained realization."""
        # set up fields
        if delta_fix is None:
            delta = self.gen_delta()
        else:
            delta = delta_fix
        t = np.copy(delta)
        tlm = hp.map2alm(t*self.weights_map, self.params.lmax, iter=0)
        # get cooling schedule
        lamb_list = self.cs.lamb_list
        eps_list = self.cs.eps_list
        for i in range(len(lamb_list)):
            flm, tlm = self.solve_flm(tlm, delta, lamb_list[i], eps_list[i])
        return flm

    # helper functions
    def check_nside(self, m):
        assert hp.get_nside(m) == self.params.nside, \
            "Wrong resolution, NSIDE should be {}.".format(self.params.nside)


def field_trnsfrm(delta, tlm, lamb, sl, tau, noise_cov, good_pix, bad_pix,
                  weights_map, pix_area, lmax, nside):
    """Transforms fluctuation field to pixel space and back given constrained
    realization, cooling schedule parameters and random field delta"""
    flm = (sl*tlm)/(sl+lamb*pix_area*tau)
    f = hp.alm2map(flm, nside, lmax, verbose=False)
    t = np.zeros(hp.nside2npix(nside))
    t[bad_pix] = lamb*tau*delta[bad_pix]+f[bad_pix]
    t[good_pix] = \
        (lamb*tau*noise_cov[good_pix]*delta[good_pix] +
         (noise_cov[good_pix] - tau)*f[good_pix])/(noise_cov[good_pix] +
                                                   (lamb-1)*tau)
    tlm = hp.map2alm(t*weights_map, lmax, iter=0)
    return flm, tlm
