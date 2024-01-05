#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 14:16:53 2022

@author: imchugh
"""

import os
import linecache
import numpy as np
# cartopy replaces basemap for plotting maps
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.geodesic as cgeo
# matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import from_levels_and_colors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import BoundaryNorm
from matplotlib_scalebar.scalebar import ScaleBar
import matplotlib.ticker as mticker
import matplotlib.colors as colors
import matplotlib as mpl
# pandas
import pandas as pd
# scipy
from scipy.signal import savgol_filter

import pdb

#------------------------------------------------------------------------------
# CLASSES #
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
class contour_mapper():

    def __init__(self, dem_asc_path):

        self.data = np.loadtxt(dem_asc_path, skiprows=6)
        self.data_header = [linecache.getline(dem_asc_path, i)
                            for i in range(1,7)]
        self.data_dict = {x.split(' ')[0]: float(x.split(' ')[1].strip())
                          for x in self.data_header}

    #--------------------------------------------------------------------------
    def get_bounds(self):

        return {'southernmost_lat': self.data_dict['yllcorner'],
                'northernmost_lat': self.data_dict['yllcorner'] +
                                    self.data_dict['cellsize'] *
                                    self.data_dict['nrows'],
                'westernmost_lon': self.data_dict['xllcorner'],
                'easternmost_lon': self.data_dict['xllcorner'] +
                                   self.data_dict['cellsize'] *
                                   self.data_dict['ncols']}
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def get_centre(self):

        d = self.get_bounds()
        return {'centre_lat': (d['southernmost_lat'] +
                               d['northernmost_lat']) / 2,
                'centre_lon': (d['easternmost_lon'] +
                               d['westernmost_lon']) / 2}
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def get_coords_as_xy(self):

        d = self.get_bounds()
        x = np.linspace(
            d['westernmost_lon'], d['easternmost_lon'], self.data.shape[1]
            )
        y = np.linspace(
            d['southernmost_lat'], d['northernmost_lat'], self.data.shape[0]
            )
        if d['southernmost_lat'] < 0:
            y = y [::-1]
        return {'x': x, 'y': y}
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def get_coord_mesh(self):

        d = self.get_coords_as_xy()
        lons, lats = np.meshgrid(d['x'], d['y'])
        return {'longitudes': lons, 'latitude': lats}
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def get_elevation_stats(self):

        return {'min': self.data.min(),
                'max': self.data.max(),
                'mean': self.data.mean()}
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def plot_map(self, colormap='Greys', marker_lat=None, marker_lon=None,
                 cmap_lower=0, cmap_upper=1, reverse_cmap=False,
                 contour_levels=None, smooth_factor=None, save_path=None):

        # Basic checks
        if not (0 >= cmap_lower <=1) and (0 >= cmap_upper <=1):
            raise ValueError('cmap bounds must be between 0 and 1')
        if not cmap_lower < cmap_upper:
            raise ValueError('cmap_lower must be less than cmap_upper')

        # Create plot and set projection and centre
        centre_dict = self.get_centre()
        fig, ax = plt.subplots(
            figsize=(20,10),
            subplot_kw=dict(projection=ccrs.Mercator(
                central_longitude=centre_dict['centre_lon'])
                )
            )

        # Set plot bounds (lower left lon, upper right lon,
        #                  lower left lat, upper right lat)
        bounds_dict = self.get_bounds()
        extent = [bounds_dict['westernmost_lon'],
                  bounds_dict['easternmost_lon'],
                  bounds_dict['southernmost_lat'],
                  bounds_dict['northernmost_lat']]
        ax.set_extent(extent)

        # Basic plot formatting
        fig.patch.set_facecolor(color=(30/255, 30/255, 30/255))
        cmap = plt.get_cmap(colormap)
        if reverse_cmap:
            new_cmap = _truncate_colormap(cmap, minval=cmap_upper,
                                          maxval=cmap_lower)
        else:
            new_cmap = _truncate_colormap(cmap, minval=cmap_lower,
                                          maxval=cmap_upper)
        stats_dict = self.get_elevation_stats()
        ax.xaxis.tick_bottom()
        ax.yaxis.tick_left()
        lon_formatter = LongitudeFormatter(dateline_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        ax.grid(linewidth=1, color='black', alpha=0.5, linestyle='--')

        # Create longitude / latitude grid onto which to project elevation data
        coords_dict = self.get_coords_as_xy()
        lons, lats = np.meshgrid(coords_dict['x'], coords_dict['y'])

        # Create elevation colour layer
        colormesh = ax.pcolormesh(lons, lats, self.data, vmin=stats_dict['min'],
                                  vmax=stats_dict['max'], cmap=new_cmap,
                                  transform=ccrs.PlateCarree(), zorder=0)
        
        # Smooth the data if required
        if smooth_factor:
            data = (
                pd.DataFrame(self.data)
                .apply(savgol_filter, window_length=smooth_factor, polyorder=3)
                .to_numpy()
                )
        else:
            data = self.data

        # Create and label contours
        if not contour_levels:
            contour_levels = _get_contour_steps(the_max=stats_dict['max'],
                                                the_min=stats_dict['min'])
        contour = ax.contour(lons, lats, data, levels=contour_levels,
                             colors='k', transform=ccrs.PlateCarree(),
                             zorder=7, alpha=0.5)
        ax.clabel(contour, inline=True, fmt='%1.0f', fontsize=12, colors='k')

        # Create colorbar
        cax, kw = mpl.colorbar.make_axes(ax, location='left',
                                         pad=0.05, shrink=0.7)
        cb = fig.colorbar(colormesh, cax=cax, **kw) # extend='both'
        cb.set_label('Elevation (m)', labelpad=5, size=20, color='White')
        cax.tick_params(labelsize=20, labelcolor='White', color='White')

        # Add features
        ax.add_feature(cfeature.RIVERS)

        # Create marker
        if not marker_lat: marker_lat = centre_dict['centre_lat']
        if not marker_lon: marker_lon = centre_dict['centre_lon']
        point_x, point_y = (marker_lon, marker_lat)
        ax.plot(point_x, point_y, marker='s', color='black', markersize=10,
                transform=ccrs.PlateCarree(), zorder=7)

        # Create marker text
        text_x = (
            centre_dict['centre_lon'] +
            (bounds_dict['easternmost_lon'] - marker_lon) / 20
            )
        text_y = marker_lat
        ax.text(text_x, text_y, 'Tower', fontsize=18,
                transform=ccrs.PlateCarree(), zorder=7, va='center')

        # Create north arrow
        arrow_lon = (
            bounds_dict['easternmost_lon'] -
            abs(bounds_dict['westernmost_lon'] -
                bounds_dict['easternmost_lon']) / 12
            )
        arrow_lat = (
            bounds_dict['northernmost_lat'] -
            abs(bounds_dict['northernmost_lat'] -
                bounds_dict['southernmost_lat']) / 15
            )
        arrow_x, arrow_y = (arrow_lon, arrow_lat)
        text_x, text_y = (0, -120)
        transform = ccrs.PlateCarree()._as_mpl_transform(ax)
        ax.annotate('N', xy=(arrow_x, arrow_y), xycoords=transform,
            xytext = (text_x, text_y), textcoords='offset points',
            color='black', fontsize=24, va='center', ha='center',
            arrowprops=dict(shrink=0.05, width=2.0, color='black'), zorder=7)

        # Create scale bar
        pt_a = (0.1, 0.1)
        pt_a_disp = ax.transAxes.transform(pt_a)
        pt_a_data = ax.transData.inverted().transform(pt_a_disp)
        pt_b = (0.1, 0.2)
        pt_b_disp = ax.transAxes.transform(pt_b)
        pt_b_data = ax.transData.inverted().transform(pt_b_disp)
        dx = (
            round((pt_b_data[1] - pt_a_data[1]) /
                  (pt_b_disp[1] - pt_a_disp[1]), 2)
            )

        d = dict(family='serif', size='xx-large')
        scalebar = ScaleBar(dx=1, units="m", location='lower left', 
                            box_alpha=0.2, pad=0.5, border_pad=4, frameon=True,
                            sep=5, scale_loc='top', font_properties=d)
        ax.add_artist(scalebar)

        # Output the image as a file
        if save_path:
            plt.savefig(save_path, bbox_inches='tight')
    #--------------------------------------------------------------------------

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# FUNCTIONS #
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def _truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n))
        )
    return new_cmap
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def _get_contour_steps(the_max, the_min, n_steps=10):

    allowed_range = np.array([5, 10, 20, 25, 30, 40, 50, 60, 75, 80, 100])
    true_range = the_max - the_min
    true_step = true_range / n_steps
    allowed_step = allowed_range[(abs(allowed_range - true_step)).argmin()]
    start = int(the_min / allowed_step) * allowed_step
    end = np.ceil(the_max / allowed_step) * allowed_step
    the_range = (
        np.linspace(start=start, stop=end,
                    num=int((end-start)/allowed_step) + 1)
        )
    return the_range
#------------------------------------------------------------------------------