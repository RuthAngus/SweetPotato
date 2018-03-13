"""
DETECTION EFFICIENCY FUNCTIONS
==============================
Functions for calculating the detection efficiency.
"""

import numpy as np


def get_duration(period, aor, e):
    """
    Equation (1) from Burke et al. This estimates the transit
    duration in the same units as the input period. There is a
    typo in the paper (24/4 = 6 != 4).

    :param period: the period in any units of your choosing
    :param aor:    the dimensionless semi-major axis (scaled
                   by the stellar radius)
    :param e:      the eccentricity of the orbit

    """
    return 0.25 * period * np.sqrt(1 - e**2) / aor

def get_a(period, mstar, Go4pi=2945.4625385377644/(4*np.pi*np.pi)):
    """
    Compute the semi-major axis of an orbit in Solar radii.

    :param period: the period in days
    :param mstar:  the stellar mass in Solar masses

    """
    return (Go4pi*period*period*mstar) ** (1./3)


def get_delta(k, c=1.0874, s=1.0187):
    """
    Estimate the approximate expected transit depth as a function
    of radius ratio. There might be a typo here. In the paper it
    uses c + s*k but in the public code, it is c - s*k:
    https://github.com/christopherburke/KeplerPORTs

    :param k: the dimensionless radius ratio between the planet and
              the star

    """
    delta_max = k*k * (c + s*k)
    return 0.84 * delta_max


def get_mes(star, period, rp, tau, cdpp_cols, cdpp_vals, re=0.009171):
    """
    Estimate the multiple event statistic value for a transit.

    :param star:   a pandas row giving the stellar properties
    :param period: the period in days
    :param rp:     the planet radius in Earth radii
    :param tau:    the transit duration in hours

    """
    # Interpolate to the correct CDPP for the duration.
    cdpp = np.array(star[cdpp_cols], dtype=float)
    sigma = np.interp(tau, cdpp_vals, cdpp)

    # Compute the radius ratio and estimate the S/N.
    k = rp * re / star.radius
    snr = get_delta(k) * 1e6 / sigma

    # Scale by the estimated number of transits.
    ntrn = star.dataspan * star.dutycycle / period
    return snr * np.sqrt(ntrn)


def get_pdet(star, aor, period, rp, e, cdpp_cols, cdpp_vals, pgam,
             mesthres_cols, mesthres_vals):
    """
    Equation (5) from Burke et al. Estimate the detection efficiency
    for a transit.

    :param star:   a pandas row giving the stellar properties
    :param aor:    the dimensionless semi-major axis (scaled
                   by the stellar radius)
    :param period: the period in days
    :param rp:     the planet radius in Earth radii
    :param e:      the orbital eccentricity

    """
    tau = get_duration(period, aor, e) * 24.
    mes = get_mes(star, period, rp, tau, cdpp_cols, cdpp_vals)
    mest = np.interp(tau, mesthres_vals,
                     np.array(star[mesthres_cols],
                              dtype=float))
    x = mes - 4.1 - (mest - 7.1)
    return pgam.cdf(x)


def get_pwin(star, period):
    """
    Equation (6) from Burke et al. Estimates the window function
    using a binomial distribution.

    :param star:   a pandas row giving the stellar properties
    :param period: the period in days

    """
    M = star.dataspan / period
    f = star.dutycycle
    omf = 1.0 - f
    pw = 1 - omf**M - M*f*omf**(M-1) - 0.5*M*(M-1)*f*f*omf**(M-2)
    msk = (pw >= 0.0) * (M >= 2.0)
    return pw * msk


def get_pgeom(aor, e):
    """
    The geometric transit probability.

    See e.g. Kipping (2014) for the eccentricity factor
    http://arxiv.org/abs/1408.1393

    :param aor: the dimensionless semi-major axis (scaled
                by the stellar radius)
    :param e:   the orbital eccentricity

    """
    return 1. / (aor * (1 - e*e)) * (aor > 1.0)
