#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 10:12:32 2019

@author: ian
"""

import numpy as np

R1 = 6371.0088

def circ_local(latitude):

    scaling = np.cos(np.radians(latitude))
    R_local = scaling * R1
    return 2 * np.pi * R_local

def earth_circ():

    return 2 *np.pi * R1

def get_arc_length_from_lat(latitude, degrees=15):

    return circ_local(latitude) * degrees / 360

def get_arc_deg_from_lat(latitude, length=100):

    return length / circ_local(latitude) * 360.0

def convert_hms_to_decimal(degrees, minutes, seconds):

    if degrees >= 360: raise RuntimeError('Degrees must be less than 360!')
    if minutes >= 60: raise RuntimeError('Minutes must be less than 60!')
    if seconds >= 60: raise RuntimeError('Seconds must be less than 60!')
    return degrees + (minutes * 60 + seconds) / 3600.0

def get_grid_box_approx(lat, lon, side_len=2.5):

    arc_degrees_along_lat = get_arc_deg_from_lat(latitude=lat, length=side_len)
    arc_degrees_along_lon = side_len / earth_circ() * 360
    return {
        'north_boundary': lat + arc_degrees_along_lon,
        'south_boundary': lat - arc_degrees_along_lon,
        'east_boundary': lon + arc_degrees_along_lat,
        'west_boundary': lon - arc_degrees_along_lat
        }
