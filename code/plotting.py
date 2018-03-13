"""
PLOTTING
========
Functions for making plots.
"""

import numpy as np
import matplotlib.pyplot as pl


def plot_obs_and_get_bin_edges(kois, var_names, nbins, ranges, log=False):
    """
    Make histograms of the observed planet properties.
    Plot the histograms and return the edges of the bins.
    params:
    --------
    var_names: (list)
    A list of strings of variable names. "period", "radius", e.g.
    nbins: (list)
    A list of the number of bins for each variable.
    ranges: (list)
    A list of tuples of variable ranges.
    log: (bool)
    Whether to make bins in linear or log space.

    returns:
    --------
    bins: (list)
    A list of lists of bin edges. One list for each variable.
    """
    nvar = len(var_names)

    # Define the bin edges.
    bins = []

    # Logarithmic.
    if log:
        for i in range(nvar):
            bins.append(np.exp(np.linspace(np.log(ranges[i][0]),
                                           np.log(ranges[i][1]), nbins[i])))

    # Linear
    else:
        for i in range(nvar):
            bins.append(np.linspace(ranges[i][0], ranges[i][1], nbins[i]))

    for i in range(nvar):
        print(i)
        hist, _ = np.histogram(kois["{}".format(var_names[i])], bins[i])
        pl.step(bins[i][:-1] + .5 * np.diff(bins[i]), hist)
        pl.xlabel("{}".format(var_names[i]))
        pl.ylabel("# planets")
        pl.show()
    return bins


def plot_marginal_completeness_2d(x_grid, y_grid, marg_comp, xlabel, ylabel,
                                  contours=False):
    """
    Make a 2d histogram plot of the completeness as a function of x and y.
    """
    pl.pcolor(x_grid, y_grid, marg_comp, cmap="BuGn")

    if contours:
        c = pl.contour(x_grid, y_grid, marg_comp, #/ len(stlr),
                       colors="k", alpha=0.8)
        pl.clabel(c, fontsize=12, inline=1, fmt="%.3f")

    pl.title("mean pipeline detection efficiency")
    pl.xlabel("{}".format(xlabel))
    pl.ylabel("{}".format(ylabel));
    pl.colorbar()
    pl.show()


def plot_det_eff_2d(grids, comp, label):
    """
    Plot the 2d histogram of the detection efficiency over 3 dimensions.
    """
    period_grid, rp_grid, c_grid = grids

    # period = 9, radius = 11, c = 10
    comp_moc = np.sum(comp, axis=2)  # Marginalize over the stellar parameter.
    comp_mor = np.sum(comp, axis=1)  # Marginalize over radius.
    comp_mop = np.sum(comp, axis=0)  # Marginalize over period.

    # Plot the three combinations
    plot_marginal_completeness_2d(period_grid[:, :, 0], rp_grid[:, :, 0],
                                  comp_moc, "Period", "$R_p \, [R_\oplus]$")
    plot_marginal_completeness_2d(period_grid[:, 0, :], c_grid[:, 0, :],
                                  comp_mor, "Period", label)
    plot_marginal_completeness_2d(rp_grid[0, :, :], c_grid[0, :, :],
                                  comp_mop, "$R_p \, [R_\oplus]$", label)


def plot_marginal_completeness_1d(x_grid, marg_comp, xlabel, ylabel)
    """
    Make a 1d histogram of the completeness as a function of x.
    """
    pl.step(x_grid, marg_comp)
    pl.xlabel("{}".format(xlabel))
    pl.ylabel(ylabel);
    pl.show()


def plot_det_eff_1d(grids, comp, var_names):
    """
    Plot the 1d histogram of the detection efficiency over 2 dimensions.
    """
    period_grid, rp_grid, c_grid = grids

    comp_c = np.sum(comp, axis=(0, 1))  # Marginalize over period and radius
    comp_p = np.sum(comp, axis=(1, 2))  # Marginalize over radius and stellar param
    comp_r = np.sum(comp, axis=(0, 2))  # Marginalize over period and stellar param

    # Plot the three combinations
    plot_marginal_completeness_1d(period_grid[:, 0, 0], comp_p, var_names[0],
                                  "Completeness")
    plot_marginal_completeness_1d(rp_grid[0, :, 0], comp_r, var_names[1],
                                  "Completeness")
    plot_marginal_completeness_1d(c_grid[0, 0, :], comp_c, var_names[2],
                                  "Completeness")


def plot_occ_rates(kois, grids, comp, var_names, bins):
    """
    Plot the observed distribution divided by the detection efficiency over 3
    dimensions.
    """
    period_grid, rp_grid, c_grid = grids

    comp_c = np.sum(comp, axis=(0, 1))  # Marginalize over period and radius
    comp_p = np.sum(comp, axis=(1, 2))  # Marginalize over radius and stellar param
    comp_r = np.sum(comp, axis=(0, 2))  # Marginalize over period and stellar param

    hist_p, _ = np.histogram(kois["{}".format(var_names[0])], bins[0])
    plot_marginal_completeness_1d(period_grid[:-1, 0, 0], hist_p/comp_p[:-1],
                                  var_names[0], "Occurrence rate")

    hist_r, _ = np.histogram(kois["{}".format(var_names[1])], bins[1])
    plot_marginal_completeness_1d(rp_grid[0, :-1, 0], hist_r/comp_r[:-1],
                                  var_names[1], "Occurrence rate")

    hist_c, _ = np.histogram(kois["{}".format(var_names[2])], bins[2])
    plot_marginal_completeness_1d(c_grid[0, 0, :-1], hist_c/comp_c[:-1],
                                  var_names[2], "Occurrence rate")
