#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 14:20:14 2024

@author: imchugh
"""

import datetime as dt
import pathlib
import warnings

import matplotlib.pyplot as plt
# import numpy as np
from cartopy import crs as ccrs #, feature as cfeature
# import cartopy.io.shapereader as shpreader
from matplotlib.image import imread

import sparql_site_details as sd

#  Suppress warnings issued by Cartopy when downloading data files
warnings.filterwarnings('ignore')
###############################################################################

###############################################################################
### CONSTANTS ###
###############################################################################
format_dict = {
    'AliceSpringsMulga': {'ha': 'right', 'xytext': (-10, 0), 'arrow': None},
    'AlpinePeat': {'ha': 'left', 'xytext': (10, -3), 'arrow': None},
    'Boyagin': {'ha': 'left', 'xytext': (25, 0), 'arrow': {'arrowstyle': '-'}},
    'Calperum': {'ha': 'left', 'xytext': (10, 0), 'arrow': None},
    'CowBay': {'ha': 'left', 'xytext': (10, 3), 'arrow': None},
    'CumberlandPlain': {'ha': 'left', 'xytext': (10, -3), 'arrow': None},
    'DalyUncleared': {'ha': 'left', 'xytext': (10, 0), 'arrow': None},
    'DryRiver': {'ha': 'left', 'xytext': (10, -3), 'arrow': None},
    'Fletcherview': {'ha': 'left', 'xytext': (10, 0), 'arrow': None},
    'Gingin': {'ha': 'left', 'xytext': (10, 0), 'arrow': None},
    'GreatWesternWoodlands': {'ha': 'left', 'xytext': (10, 0), 'arrow': None},
    'HowardSprings': {'ha': 'right', 'xytext': (-10, 0), 'arrow': None},
    'Litchfield': {'ha': 'left', 'xytext': (10, 0), 'arrow': None},
    'Longreach': {'ha': 'left', 'xytext': (10, 0), 'arrow': None},
    'MyallValeA': {'ha': 'left', 'xytext': (10, 0), 'arrow': None},
    'Ridgefield': {'ha': 'left', 'xytext': (10, -20), 'arrow': {'arrowstyle': '-'}},
    'RobsonCreek': {'ha': 'left', 'xytext': (10, -3), 'arrow': None},
    'Samford': {'ha': 'left', 'xytext': (10, 0), 'arrow': None},
    'SilverPlain': {'ha': 'left', 'xytext': (10, 5), 'arrow': None},
    'SturtPlains': {'ha': 'left', 'xytext': (10, 0), 'arrow': None},
    'TiTreeEast': {'ha': 'left', 'xytext': (10, 0), 'arrow': None},
    'Tumbarumba': {'ha': 'left', 'xytext': (10, 0), 'arrow': None},
    'Warra': {'ha': 'left', 'xytext': (10, -10), 'arrow': None},
    'Wellington': {'ha': 'left', 'xytext': (10, 3), 'arrow': None},
    'Whroo': {'ha': 'right', 'xytext': (-8, 3), 'arrow': None},
    'WombatStateForest': {'ha': 'right', 'xytext': (-8, -3), 'arrow': None},
    'Yanco': {'ha': 'left', 'xytext': (5, 10), 'arrow': None}
    }

legend_dict = {
    'green': '<1 day',
    'yellow': '$1\leq\/days<3$',
    'orange': '$3\leq\/days<5$',
    'magenta': '$5\leq\/days<7$',
    'red': '$\geq\/7 days$',
    'grey': 'N/A'
    }
###############################################################################

###############################################################################
### MAIN FUNCTION ###
###############################################################################

def main():

    # Set projection and extent
    projPC = ccrs.PlateCarree()
    lonW = 110
    lonE = 158
    latS = -45
    latN = -8
    res = '50m'
    
    # Set background and border shape file locations
    bkgnd_file = pathlib.Path(
        '/home/unimelb.edu.au/imchugh/Dropbox/Work/Data/Background_images/'
        'NE1_50M_SR_W.tif'
        )
    
    # This is busted ATM
    # borders_file = pathlib.Path(
    #     '/home/unimelb.edu.au/imchugh/Dropbox/Work/Data/Layers/gadm41_AUS_0.shp'
    #     )
    # adm1_shapes = list(shpreader.Reader(borders_file).geometries())
    
    deets = sd.site_details()
    
    # Plot
    fig = plt.figure(figsize=(11, 8.5))
    ax = plt.subplot(1, 1, 1, projection=projPC)
    ax.set_title(
        'Network connection status: '
        f'{dt.datetime.now().strftime("%Y-%m-%d %H:%M")} AEST',
        fontsize=18,
        pad=20
        )
    gl = ax.gridlines(
        draw_labels=True, linewidth=2, color='gray', alpha=0, linestyle='--'
        )
    ax.coastlines(resolution=res, color='black')
    ax.imshow(imread(bkgnd_file), origin='upper', transform=projPC, 
              extent=[-180, 180, -90, 90])
    ax.set_extent([lonW, lonE, latS, latN], crs=projPC)
    
    for item in legend_dict.items():
        ax.plot(
            0.5, -0.5, transform=ax.transAxes, marker='o', ms=5, color=item[0], 
            linestyle='None', label=item[1],
            )
    
    for site in deets.get_operational_sites().index:
        site_deets = deets.get_single_site_details(site=site)    
        ax.plot(
            site_deets.longitude, 
            site_deets.latitude, 
            transform=ccrs.Geodetic(), marker='o',
            color='grey', ms=10, alpha=0.6
            )
        ax.annotate(
            site, xy=(site_deets.longitude, site_deets.latitude), 
            xytext=format_dict[site]['xytext'], 
            textcoords='offset points', ha=format_dict[site]['ha'], 
            va='center', 
            arrowprops=format_dict[site]['arrow'], 
            transform=ccrs.Geodetic()
            )
    
    ax.legend(
        title='Time since last data', ncols=2, loc=(0.06, 0.048), markerscale=1.5
        )
###############################################################################

if __name__=='__main__':
    
    main()