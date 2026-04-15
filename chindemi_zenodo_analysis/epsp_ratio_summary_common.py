import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize, stats
from scipy.interpolate import griddata
from matplotlib import cm, colors
from itertools import product
from matplotlib.transforms import Affine2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# Common functions
def plot_ratio(df0, label0, df1=None, label1=None,
               style='match', xvar='frequency_train',
               offset=True, marker=None, ax=None):
    """Plot EPSP ratio in dataset df0 as a function of xvar (i.e. stimulation frequency).
    Two different datasets can be compared by specifying an additional df1. If df1 is provided, the
    style parameter allows to select the type of comparison ('match' or 'diverge'). For the moment
    the only difference bethween the two styles is the color used to disply elements in df1.
    Furthermore, the means of df0 and df1 will be compared using the Welch's t-test and the results
    printed on the screen.
    
    Parameters
    ----------
    df0 : pandas.DataFrame
        Main dataset to use for the plot
    label0 : str
        Name of the main dataset, used in the legend
    df1 : pandas.DataFrame
        Optional extra dataset to use for the plot
    label1 : str
        Name of the extra dataset, used in the legend
    style : str ('match' or 'diverge')
        Plot style of the comparison.
    xvar : str
        Column to use as x variable in the plot
    offset : bool
        Apply a small offset to the x coordinate to avoid the overlap of the errorbars
    marker : char
        Marker for EPSP ratio mean.
    ax : matplotlib.axes.Axes
        Plot on existing axes

    Returns
    -------
    fig : matplotlib.Figure
    ax : matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    if offset:
        trans1 = Affine2D().translate(-0.7, 0.0) + ax.transData
        trans2 = Affine2D().translate(+0.7, 0.0) + ax.transData
    else:
        trans1 = None
        trans2 = None

    if marker is None:
        fmt = 'o-'
    else:
        fmt = marker + '-'

    # Plot main
    df0_s = df0.sort_values(xvar)
    ax.errorbar(df0_s[xvar], df0_s.mean_epsp_ratio, yerr=df0_s.sem_epsp_ratio,
                fmt=fmt, color='C0', transform=trans1, label=label0)

    # Plot secondary
    if df1 is not None:
        df1_s = df1.sort_values(xvar)
    
        if style == 'match':
            ax.errorbar(df1_s[xvar], df1_s.mean_epsp_ratio, yerr=df1_s.sem_epsp_ratio,
                        fmt=fmt, color='C1', transform=trans2, label=label1)
        elif style == 'diverge':
            ax.errorbar(df1_s[xvar], df1_s.mean_epsp_ratio, yerr=df1_s.sem_epsp_ratio,
                        fmt=fmt, color='lightgrey', transform=trans2, label=label1)
        else:
            raise NotImplemented('Unknown plot style')
            
        # Welch’s t-test
        result = pd.merge(df0_s, df1_s, on=xvar)
        wt_stat, wt_p = stats.ttest_ind_from_stats(
                result.mean_epsp_ratio_x, result.std_epsp_ratio_x, result.sample_size_x,
                result.mean_epsp_ratio_y, result.std_epsp_ratio_y, result.sample_size_y,
                equal_var=False)
        for x in zip(result[xvar], wt_stat, wt_p):
            print("t-test results for %s = %f: s = %.3f, p = %.3f" % (xvar, *x))

    # Finalize plot
    if xvar == 'frequency_train':
        ax.set_xlabel('Frequency (Hz)')
    elif xvar == 'dt_train':
        ax.axvline(0, color='k', lw=.5, ls='dashed')
        ax.set_xlabel(r'$\Delta t$ (ms)')
    ax.set_ylabel('EPSP ratio')
    ax.axhline(1, color='k', lw=.5, ls='dashed')
    ax.legend(loc=0)
    
    return fig, ax


def plot_map(df, vmin=.7, vmax=1.45, mapsize=None):
    """Plot EPSP ratio for each pair of stimulation frequency and pre-post timing in df.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Main dataset to use for the plot
    vmin : float
        Minimum EPSP ratio to display in the colorbar
    vmax : float
        Maximum EPSP ratio to display in the colorbar

    Returns
    -------
    fig : matplotlib.Figure
    """
    grid_x, grid_y = np.mgrid[df.dt_train.min():df.dt_train.max():1000j,
                             df.frequency_train.min():df.frequency_train.max():1000j]
    grid_z = griddata(df[['dt_train', 'frequency_train']], df.mean_epsp_ratio,
                      (grid_x, grid_y), method='cubic')
    fig, ax = plt.subplots(figsize=mapsize)
    ax.plot(df.dt_train, df.frequency_train, 'k.', ms=1)
    assert vmin < np.nanmin(grid_z)
    assert vmax > np.nanmax(grid_z)
    divnorm = colors.TwoSlopeNorm(vmin=vmin, vcenter=1, vmax=vmax)
    extent = (df.dt_train.min(),df.dt_train.max(),df.frequency_train.min(),df.frequency_train.max())
    pos = ax.imshow(grid_z.T, cmap=cm.turbo, origin='lower',
                    norm=divnorm, aspect='auto', extent=extent)
    fig.colorbar(pos, extend='both', ax=ax, label='EPSP ratio')
    ax.set_ylabel('Frequency (Hz)')
    ax.set_xlabel(r'$\Delta t$ (ms)')
    return fig


def load_insilico(datapath, epsp_cut=0.01):
    """Load in silico results and compute aggregated stats.
    The function filters out outliers, based on initial EPSP size.
    Eliminated outliers are printed on the screen.
    
    Parameters
    ----------
    datapath : str
        Path of CSV results file
    epsp_cut : float
        Minimum initial EPSP amplitude allowed
    
    Return
    ------
    data_clean : pandas.DataFrame
        Results of individual experiments, after outliers removal
    data_agg : pandas.DataFrame
        Aggregated results over experimental conditions (stimulation frequency and timing)
    """
    data = pd.read_csv(datapath)

    # Find EPSP amplitude column
    if "epsp_before_amplitude" in data.columns:
        # Slope dataset, EPSP amplitude was explicitly added
        epsp_before_col = "epsp_before_amplitude"
    else:
        # Amplitude dataset
        epsp_before_col = "epsp_before"

    # Detect outliers (unrealistically weak connections)
    outliers = data[data[epsp_before_col] < epsp_cut]
    if outliers.shape[0] > 0:
        print('Found outliers in %s:' % datapath)
        display(outliers[['pregid', 'postgid', 'frequency_train',
                          'dt_train', 'epsp_ratio', epsp_before_col]])
    else:
        print('No outliers in %s' % datapath)
    data_clean = data.drop(outliers.index)

    # Compute aggregated stats
    data_agg = (data_clean
                .groupby(['frequency_train', 'dt_train'])['epsp_ratio']
                .agg(['mean', 'std', 'sem', 'count'])
                .add_suffix('_epsp_ratio')
                .reset_index())
    data_agg.dropna(inplace=True)
    data_agg.rename(columns={'count_epsp_ratio': 'sample_size'}, inplace=True)

    print("Data sample:")
    display(data_clean.sample(5, random_state=1234)[['pregid', 'postgid', 'frequency_train',
                                                     'dt_train', 'epsp_ratio', epsp_before_col]])
    
    print("Aggregated data:")
    display(data_agg)

    return data_clean, outliers, data_agg


def powerlaw(x, a, b):
    return a*x**b


def fit_powerlaw(x, y):
    """Fit the powerlaw a*x^b.
    The fit is perfomed in logaritmic space. Propagation of errors is used to estimate uncertainty
    on the *a* parameter.
    
    See also https://scipy-cookbook.readthedocs.io/items/FittingData.html
    """
    logx = np.log10(x)
    logy = np.log10(y)
    
    popt, pcov = optimize.curve_fit(lambda x, c, b: c + b*x, logx, logy)
    perr = np.sqrt(np.diag(pcov))
    
    a = 10**popt[0]
    b = popt[1]
    
    a_error = perr[0] * a
    b_error = perr[1]

    return (a, b), (a_error, b_error)
