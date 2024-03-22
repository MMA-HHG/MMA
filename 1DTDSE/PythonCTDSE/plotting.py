"""
Plotting
========

This module contains methods for plotting the PythonCTDSE results conveniently.

"""
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
 
def plot_colormap(
        y_axis, 
        x_axis, 
        z, 
        plot_scale = "log", 
        z_min = 1e-14,
        z_max = 1,
        figsize = (12, 3),
        y_label = "",
        x_label = "",
        cmap = "jet"
    ):
    """
    Plots the colormap of arbitrary matrix for arbitrary distribution.
    
    Parameters:
    -----------
    y_axis (numpy.ndarray):
        Values on the y-axis on the colormap.
    x_axis (numpy.ndarray):
        Values on the y-axis on the colormap.
    z (numpy.ndarray):
        Matrix for plotting.
    plot_scale: {'log', 'linear'}, optional, default: 'log'
        Plot with logarithmic or linear scale.
    z_min (float):
        Minimum plotting value on the colormesh.
    z_max (float):
        Maximum threshold for the colorbar
    figsize (tuple), optional, default: (12, 3)
        Set figure size in inches - (width, height)
    y_label (str), optional, default: ""
        Label of the y-axis
    x_label (str), optional, default: ""
        Label of the x-axis
    cmap (str), optional, default: "jet"
        Colormap name from the Matplotlib colormap library

    Returns:
    --------
    None
    """
    fig, ax = plt.subplots()
    
    if plot_scale == 'log':
        col = LogNorm(vmin = z_min, vmax = z_max)
    elif plot_scale == 'linear':
        col = Normalize(vmin = z_min, vmax = z_max)
    else:
        raise ValueError("Unknown scale '" + plot_scale +"'. "
                         "Available options are 'log' and 'linear'")
        
    c = ax.pcolormesh(x_axis, y_axis, z,
        cmap = cmap,
        shading = 'gouraud',
        norm = col
    )

    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)

    fig.colorbar(c, ax=ax)
    fig.dpi = 300
    fig.set_size_inches(figsize)

    plt.show()

def plot(
        *args,
        plot_scale = "linear", 
        figsize = (5, 3),
        y_label = "",
        x_label = "",
        **kwargs
        ):
    """
    Simple plot method.

    Parameters:
    -----------
    *args
        Arguments for the plot command, i.e. x, y or x1, y1, x2, y2, .. .
    plot_scale: {'log', 'linear'}, optional, default: 'linear'
        Plot with logarithmic or linear scale.
    figsize (tuple), optional, default: (5, 3)
        Set figure size in inches - (width, height)
    y_label (str), optional, default: ""
        Label of the y-axis
    x_label (str), optional, default: ""
        Label of the x-axis
    **kwargs
        Keyword arguments for the plot command.

    Returns:
    --------
    None
    """
    
    fig, ax = plt.subplots()

    if plot_scale == 'log':
        plt.semilogy(*args, **kwargs)
    elif plot_scale == 'linear':
        plt.plot(*args, **kwargs)
    else:
        raise ValueError("Unknown scale '" + plot_scale +"'. "
                         "Available options are 'log' and 'linear'")
    fig.dpi = 300
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    fig.set_size_inches(figsize)
    
    plt.show()