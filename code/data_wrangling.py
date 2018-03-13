"""
DATA WRANGLING
==============
Functions for loading data and making cuts.
"""

import os
import numpy as np
import pandas as pd
import requests
from io import StringIO, BytesIO
import matplotlib.pyplot as pl


def get_catalog(name, basepath="data"):
    fn = os.path.join(basepath, "{0}.h5".format(name))
    if os.path.exists(fn):
        return pd.read_hdf(fn, name)
    if not os.path.exists(basepath):
        os.makedirs(basepath)
    print("Downloading {0}...".format(name))
    url = ("http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/"
           "nph-nstedAPI?table={0}&select=*").format(name)
    r = requests.get(url)
    if r.status_code != requests.codes.ok:
        r.raise_for_status()
    fh = BytesIO(r.content)
    df = pd.read_csv(fh)
    df.to_hdf(fn, name, format="t")
    return df


def selection(stlr, ranges):
    print("WARNING: period_rng and radius_rng should be the first and" \
          " second elements of ranges")

    # Select G and K dwarfs.
    m = (4200 <= stlr.teff) & (stlr.teff <= 6100)
    m &= stlr.radius <= 1.15

    # Only include stars with sufficient data coverage.
    m &= stlr.dataspan > 365.25*2.
    m &= stlr.dutycycle > 0.6
    m &= stlr.rrmscdpp07p5 <= 1000.

    # Only select stars with mass estimates.
    m &= np.isfinite(stlr.mass)

    base_stlr = pd.DataFrame(stlr)
    stlr = pd.DataFrame(stlr[m])

    print("Selected {0} targets after cuts".format(len(stlr)))

    kois = get_catalog("q1_q16_koi")

    # Join on the stellar list.
    kois = pd.merge(kois, stlr, on="kepid", how="inner")

    # Only select the KOIs in the relevant part of parameter space.
    m = kois.koi_pdisposition == "CANDIDATE"
    base_kois = pd.DataFrame(kois[m])
    m &= (ranges[0][0] <= kois.koi_period) & (kois.koi_period <= ranges[0][1])
    m &= np.isfinite(kois.koi_prad) & (ranges[1][0] <= kois.koi_prad) & \
        (kois.koi_prad <= ranges[1][1])

    kois = pd.DataFrame(kois[m])

    print("Selected {0} KOIs after cuts".format(len(kois)))

    yerr = np.abs(np.array(base_kois[["koi_prad_err2", "koi_prad_err1"]])).T
    pl.errorbar(base_kois.koi_period, base_kois.koi_prad, yerr=yerr, fmt=".k",
                ms=4, capsize=0, alpha=0.3)
    pl.plot(kois.koi_period, kois.koi_prad, ".k", ms=6)
    pl.fill_between(ranges[0], [ranges[1][1], ranges[1][1]],
                    [ranges[1][0], ranges[1][0]], color="g", alpha=0.2)
    pl.xlim(ranges[0] + 10 * np.array([-1, 1]))
    pl.ylim(ranges[1] + 0.5 * np.array([-1, 1]))
    pl.xlabel("period [days]")
    pl.ylabel("$R_p \, [R_\oplus]$");
    pl.show()

    return stlr, kois
