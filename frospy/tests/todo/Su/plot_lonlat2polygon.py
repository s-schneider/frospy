#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 17:29:23 2020

@author: talavera
"""

# Use with geopandas enviroment to create shapefile
# conda activate geopandas    
# conda deactivate   

import glob
import numpy                               as np
import geopandas as gpd
from geopandas import GeoSeries
from shapely.geometry import Polygon
import matplotlib.pyplot as plt

# # From Creasy et al 2019: D" anisotropy
# mdir = "/data/talavera/notes/Dst/ani_polygons/*.csv"
# ifiles = glob.glob(mdir)

# fig, ax = plt.subplots(1, 1)

# world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
# world.plot(ax=ax)


# # lonlat = sorted(lonlat, key=lambda x: x[0])
# # lonlat = list(zip(*lonlat))
# for i,f in enumerate(ifiles):
#     name = f.split(".csv")[0]
#     f = open(f, 'r') 
#     lines = f.read().splitlines()
#     lines = lines[1::]
#     lon = []
#     lat = []
#     for l in lines:
#         l = l.split()
#         lon.append(float(l[0]))
#         lat.append(float(l[1]))
#     # fig1, ax1 = plt.subplots(1, 1)
#     # print("polygon", i, len(lon), len(lat))
#     p1 =  Polygon(zip(lon, lat))
#     crs = {'init': 'epsg:4326'}
#     polygon = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[p1])       
#     polygon.geometry.plot(ax=ax,color='None', edgecolor='red')
#     # polygon.geometry.plot(ax=ax1)
#     # ax1.set_title(name)
#     polygon.to_file(filename='%s.shp'%name, driver="ESRI Shapefile")


# From Yu & GArnero, 2017:ULVZ
mdir = "/data/talavera/notes/Dst/ulvz_polygons/All.yes.grd.txt"
ifiles = glob.glob(mdir)

fig, ax = plt.subplots(1, 1)

world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
world.plot(ax=ax)


# lonlat = sorted(lonlat, key=lambda x: x[0])
# lonlat = list(zip(*lonlat))
for i,f in enumerate(ifiles):
    name = f.split(".txt")[0]
    f = open(f, 'r') 
    lines = f.read().splitlines()
    lines = lines[1::]
    lon = []
    lat = []
    for l in lines:
        l = l.split()
        if int(l[2]) == 1:
            # lon.append(float(l[0]))
            # lat.append(float(l[1]))
            ax.scatter(float(l[0]),float(l[1]), color="r", s=0.1)
    # fig1, ax1 = plt.subplots(1, 1)
    # print("polygon", i, len(lon), len(lat))
    # p1 =  Polygon(zip(lon, lat))
    # crs = {'init': 'epsg:4326'}
    # polygon = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[p1])       
    # polygon.geometry.plot(ax=ax,color='None', edgecolor='red')
    # # polygon.geometry.plot(ax=ax1)
    # # ax1.set_title(name)
    # # polygon.to_file(filename='%s.shp'%name, driver="ESRI Shapefile")
