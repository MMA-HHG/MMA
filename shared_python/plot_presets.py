import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py
import sys
import units
import mynumerics as mn

import re
import glob

import warnings

import copy


import matplotlib.pyplot as plt
from matplotlib import rcParams, rcParamsDefault


class figure_driver:
    def __init__(self):
        self.sf = []
        
        self.legend_args = []
        self.legend_kwargs = {}
        self.legend_set_in_layout = None
        
        self.right_axis_legend_args = []
        self.right_axis_legend_kwargs = {}
        
        self.savefig_args = None
        self.savefig_kwargs = {}

        self.savefigs_args = None
        self.savefigs_kwargs = []
        
        self.show_fig = True
        
        self.xlim_args = []; self.xlim_kwargs = {}        
        self.ylim_args = []; self.ylim_kwargs = {}
        self.right_ylim_args = []; self.right_ylim_kwargs = {}
        self.invert_xaxis = False
        
        # self.set_aspect_args = []; self.set_aspect_kwargs = {}  
        self.set_size_inches_args = []; self.set_size_inches_kwargs = {}  
        
        self.add_right_y_axis = False
        
        self.args = []
        self.kwargs = {'constrained_layout' : True}
        
        self.ax_minor_ticks = False
        self.ax_tick_params_kwargs = {}
        
        self.ax_right_minor_ticks = False
        self.ax_right_tick_params_kwargs = {}
        
        self.yscale = ''
        self.yscale_kwargs = {}
        
        self.set_fontsizes = ''
        
        self.annotation = None
        
        


class colorbar:
    def __init__(self):
        self.show = False
        self.show_contours = False
        self.kwargs = {}
        
class plotter:
    def __init__(self):
        self.method = plt.plot
        self.args = []
        self.kwargs = {}
        self.colorbar = colorbar()


def plot_preset(i):
    rcParams.update(rcParamsDefault)
    # param1 = copy.deepcopy(rcParams)
    if (len(i.set_fontsizes) > 0):
        if (i.set_fontsizes == 'triplet'):
            SMALL_SIZE = 17
            MEDIUM_SIZE = 20
            BIGGER_SIZE = 23
            
            plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
            plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
            plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
            plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
            plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
            plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
            plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
            
        elif (i.set_fontsizes == 'doublet'):
            SMALL_SIZE = 11.75
            MEDIUM_SIZE = 12.75
            BIGGER_SIZE = 13.75
            
            plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
            plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
            plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
            plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
            plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
            plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
            plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
            
        elif (i.set_fontsizes == 'doublet+'):
            SMALL_SIZE = 13.0
            MEDIUM_SIZE = 14.0
            BIGGER_SIZE = 15.5
            
            # SMALL_SIZE = 9.0
            # MEDIUM_SIZE = 9.0
            # BIGGER_SIZE = 30.0
            
            plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
            plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
            plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
            plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
            plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
            plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
            plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
 
        elif (i.set_fontsizes == 'largefonts'):
            SMALL_SIZE = 18.0
            MEDIUM_SIZE = 20.0
            BIGGER_SIZE = 20.5
            
            # SMALL_SIZE = 9.0
            # MEDIUM_SIZE = 9.0
            # BIGGER_SIZE = 30.0
            
            plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
            plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
            plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
            plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
            plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
            plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
            plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
            
        else: warnings.warn('fonsizes wrongly specified, default used')
    
    # param2 = copy.deepcopy(rcParams)
 
    if not(i.annotation is None):
      if not(i.annotation[0] is None):
        i.sf.append(plotter)
        i.sf[-1].method = plt.annotate  
        i.sf[-1].args = i.annotation[0]
        i.sf[-1].kwargs = i.annotation[1]
        
    Nsf = len(i.sf)
    
    fig, ax = plt.subplots(*i.args,**i.kwargs)
    
    add_legend = False
    
    right_axis_exists = False
    add_right_axis_legend = False
    
    
    for k1 in range(Nsf):
        # plot
        if (i.sf[k1].method is plt.plot):
            # choose left or right axis
            if hasattr(i.sf[k1], 'axis_location'):
              if i.sf[k1].axis_location == 'right':
                if not(right_axis_exists): ax_right = ax.twinx(); right_axis_exists = True
                
                ax_right.plot(*i.sf[k1].args,**i.sf[k1].kwargs)
                if ('label' in i.sf[k1].kwargs.keys()): add_right_axis_legend = True
                
            else:
                ax.plot(*i.sf[k1].args,**i.sf[k1].kwargs)
                if ('label' in i.sf[k1].kwargs.keys()): add_legend = True

        elif (i.sf[k1].method is plt.semilogy):
            # choose left or right axis
            if hasattr(i.sf[k1], 'axis_location'):
              if i.sf[k1].axis_location == 'right':
                if not(right_axis_exists): ax_right = ax.twinx(); right_axis_exists = True
                
                ax_right.semilogy(*i.sf[k1].args,**i.sf[k1].kwargs)
                if ('label' in i.sf[k1].kwargs.keys()): add_right_axis_legend = True
                
            else:
                ax.semilogy(*i.sf[k1].args,**i.sf[k1].kwargs)
                if ('label' in i.sf[k1].kwargs.keys()): add_legend = True
                
        elif (i.sf[k1].method is plt.errorbar):
            # choose left or right axis
            if hasattr(i.sf[k1], 'axis_location'):
              if i.sf[k1].axis_location == 'right':
                if not(right_axis_exists): ax_right = ax.twinx(); right_axis_exists = True
                
                ax_right.errorbar(*i.sf[k1].args,**i.sf[k1].kwargs)
                if ('label' in i.sf[k1].kwargs.keys()): add_right_axis_legend = True
                
            else:
                ax.errorbar(*i.sf[k1].args,**i.sf[k1].kwargs)
                if ('label' in i.sf[k1].kwargs.keys()): add_legend = True
            
        elif (i.sf[k1].method is plt.pcolor):
            map1 = ax.pcolor(*i.sf[k1].args,**i.sf[k1].kwargs)
            if i.sf[k1].colorbar.show:
                colorbar = fig.colorbar(mappable=map1, ax = ax, **i.sf[k1].colorbar.kwargs)

        elif (i.sf[k1].method is plt.pcolormesh):
            map1 = ax.pcolor(*i.sf[k1].args,**i.sf[k1].kwargs)
            if i.sf[k1].colorbar.show:
                colorbar = fig.colorbar(mappable=map1, ax = ax, **i.sf[k1].colorbar.kwargs)
                
        elif (i.sf[k1].method is plt.contourf):
            map1 = ax.contourf(*i.sf[k1].args,**i.sf[k1].kwargs)
            if i.sf[k1].colorbar.show:
                colorbar = fig.colorbar(mappable=map1, **i.sf[k1].colorbar.kwargs)
                
        elif (i.sf[k1].method is plt.contour): # if colorbar modified, has to be applied after pcolor
            map2 = ax.contour(*i.sf[k1].args,**i.sf[k1].kwargs)
            if i.sf[k1].colorbar.show_contours:
                try:
                    colorbar.add_lines(map2)
                except:
                    raise(ValueError('cannot add contours to colorbar'))
            
            # if hasattr(i.sf[k1], 'colorbar'):
            #     fig.colorbar(map1)
        elif (i.sf[k1].method is plt.hist):
            ax.hist(*i.sf[k1].args,**i.sf[k1].kwargs)
            
        elif (i.sf[k1].method is plt.annotate):
            ax.annotate(*i.sf[k1].args,**i.sf[k1].kwargs)
            
        elif (i.sf[k1].method is plt.text):
            ax.text(*i.sf[k1].args,**i.sf[k1].kwargs)
            
        elif (i.sf[k1].method is plt.vlines):
            ax.vlines(*i.sf[k1].args,**i.sf[k1].kwargs)
            
        elif (i.sf[k1].method is plt.axvline):
            ax.axvline(*i.sf[k1].args,**i.sf[k1].kwargs)
            
        elif (i.sf[k1].method is None):
            pass
                
        else:
            raise(NotImplementedError())

    if ((len(i.legend_args) > 0) or (len(i.legend_kwargs) > 0)): add_legend = True
    
    if add_legend:
        leg = ax.legend(*i.legend_args, **i.legend_kwargs)
        if not(i.legend_set_in_layout is None): leg.set_in_layout(i.legend_set_in_layout)
        
    if add_right_axis_legend:
        ax_right.legend(*i.right_axis_legend_args, **i.right_axis_legend_kwargs)
        
        
    if ((len(i.xlim_args) > 0) or (len(i.xlim_kwargs) > 0)):
        ax.set_xlim(*i.xlim_args, **i.xlim_kwargs)
        
    if i.invert_xaxis: ax.invert_xaxis()
        
        
    if ((len(i.ylim_args) > 0) or (len(i.ylim_kwargs) > 0)):
        ax.set_ylim(*i.ylim_args, **i.ylim_kwargs)
        
    if (((len(i.right_ylim_args) > 0) or (len(i.right_ylim_kwargs) > 0)) and (right_axis_exists or i.add_right_y_axis)):
        if not(right_axis_exists): ax_right = ax.twinx(); right_axis_exists = True
        ax_right.set_ylim(*i.right_ylim_args, **i.right_ylim_kwargs)
 
     
    if i.ax_minor_ticks:   
        ax.minorticks_on()
        if (len(i.ax_tick_params_kwargs) > 0): 
            ax.tick_params(**i.ax_tick_params_kwargs)
            
    if (i.ax_right_minor_ticks and right_axis_exists):   
        ax_right.minorticks_on()
        if (len(i.ax_right_tick_params_kwargs) > 0): 
            ax_right.tick_params(**i.ax_right_tick_params_kwargs)
        
    
    if (len(i.yscale) > 0):
        ax.set_yscale(i.yscale,**i.yscale_kwargs)
        
    
    if hasattr(i, 'xlabel'):
        ax.set_xlabel(i.xlabel)
    if hasattr(i, 'right_ylabel'):
        if not(right_axis_exists): ax_right = ax.twinx(); right_axis_exists = True
        ax_right.set_ylabel(i.right_ylabel)
    if hasattr(i, 'ylabel'):   
        ax.set_ylabel(i.ylabel)
    if hasattr(i, 'title'):   
        ax.set_title(i.title)
        # fig.suptitle(i.title)
        
        
    # if ((len(i.set_aspect_args) > 0) or (len(i.set_aspect_kwargs) > 0)):
    #     ax.set_aspect(*i.set_asppect_args,**i.set_asppect_args)
    
    # fig.set_size_inches(20, 10)
    if ((len(i.set_size_inches_args) > 0) or (len(i.set_size_inches_kwargs) > 0)):
        fig.set_size_inches(*i.set_size_inches_args, **i.set_size_inches_kwargs)
        
    
        
    if not(i.savefig_args is None):
        fig.savefig(*i.savefig_args, **i.savefig_kwargs)
        
    if not(i.savefigs_args is None):
        for k1 in range(len(i.savefigs_args)):
            fig.savefig(*i.savefigs_args[k1], **i.savefigs_kwargs[k1])
    
    if i.show_fig:
        plt.show()
        
    return fig
            
        