import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.spatial import cKDTree
import coordinateSystems as cs
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import glob
import interpolation
import matplotlib.path as mpth
from cartopy.io.shapereader import Reader
import json
import geojson
import rasterio
import rioxarray as rxr
import geopandas as gpd
from shapely.geometry import mapping
from metpy.units import units
import metpy.calc as mpcalc
import pandas as pd
import os


def qv_from_td(t,td,p):
    '''INPUT:
       temperature and dewpoint in Kelvin
       surface pressure in Pa

       RETURN: mixing ratio in g/kg'''
    t = t*units.K
    td = td*units.K
    p  = p*units.pascal
    rh = mpcalc.relative_humidity_from_dewpoint(t,td)
    qv = mpcalc.mixing_ratio_from_relative_humidity(p,t,rh)
    return (qv.magnitude)*1000

def centers_to_edges(xcenter,ycenter):
    edgex = np.zeros(len(xcenter)+1)
    edgey = np.zeros(len(ycenter)+1)
    edgex[0] = xcenter[0]-(xcenter[1]-xcenter[0])/2
    edgey[0] = ycenter[0]-(ycenter[1]-ycenter[0])/2
    edgex[1:] = xcenter+(xcenter[1]-xcenter[0])/2
    edgey[1:] = ycenter+(ycenter[1]-ycenter[0])/2
    return edgex,edgey

def bilinear_interp(grid1x, grid1y, grid2x, grid2y, z):
    """
    A method which interpolates a function
    z(grid1x, grid1y) of a grid (grid1x, grid1y) to another
    grid (grid2x, grid2y). Returns an array from the approximated
    function of the second grid (approximation of z(grid2x, grid2y)).

    Written by Austin Coleman, Ph.D.
    """
    # Pair flattened x and y values as coordinates
    coords_from = list(zip(grid1y.flatten(), grid1x.flatten()))
    Z = z.flatten()
    # Set up interpolation function with original grid and interp variable
    interp = interpolate.LinearNDInterpolator(coords_from, Z, fill_value=np.nan)
    # Interpolate to new grid
    interpolated_z = interp(grid2y, grid2x)
    return interpolated_z

def linear_interpolation(data, lat_source, lon_source, lat_target, lon_target):
    '''
    A function that interpolates data from a source grid onto a target grid
    all LON/LAT arrays must be 2-Dimensional (i.e., meshgrid)
    IMPORTANT: lat,lon 2-D arrays must have dimensions (lon_dim,lat_dim). 
    They must be the same between the source and target. Otherwise, you get weird columns. 
    
    '''
    # Initialize geographic coordinate system
    geo = cs.GeographicSystem()
    # Convert to Earth-Center-Earth-Fixed for common origin
    cecefx, cecefy, cecefz = geo.toECEF(lon_source, lat_source, lat_source)
    fecefx, fecefy, fecefz = geo.toECEF(lon_target, lat_target, lat_target)
    # interpolate ERA precip to WRF grid
    interp_data = bilinear_interp(cecefx, cecefy, fecefx, fecefy, data)
    return interp_data

def accadia_interpolation(lon_target,lat_target,lon_source,lat_source,data,search_radius=12500):
    '''
    Interpolates data to a new grid using a remapping technique (a simple nearest-neighbor average method)
    documented in Accadia et al. (2003). Used primarily for precipitation accumulation data.
    ----
    interp.f (interpolation.py) must be compiled with f2py before use of this function
    Originally written by Russell Manser

    INPUTS:
    data: source data to be interpolated to new grid
    lon_target,lat_target: (MxN) 2d arrays of grid to which data is interpolated
    lon_source,lat_source: 2d arrays of grid from which data is already on
    search_radius: gridspacing (plus 500m buffer) of source data grid

    OUTPUT: (MxN) reshaped array of interpolated data
    '''
    if lon_target.min()<0:
        lon_target = lon_target+360
    if lon_source.min()<0:
        lon_source = lon_source+360

    geo = cs.GeographicSystem()
    source_x,source_y,_ = geo.toECEF(lon_source,lat_source,np.zeros_like(lon_source))
    target_x,target_y,_ = geo.toECEF(lon_target,lat_target,np.zeros_like(lon_target))

    source_coords = np.array([[x,y] for x,y in zip(source_x.flatten(),source_y.flatten())])
    target_coords = np.array([[x,y] for x,y in zip(target_x.flatten(),target_y.flatten())])

    tree = cKDTree(source_coords)
    query= tree.query_ball_point(target_coords,search_radius)

    x_deltas = np.gradient(target_x)[0].flatten()
    y_deltas = np.gradient(target_y)[1].flatten()

    interp_data = interpolation.accadia_f(source_coords,data.flatten(),target_coords,x_deltas,y_deltas,query)
    interp_data.shape = target_x.shape
    return interp_data


def mask_inside_polygon_path(lon,lat,data,datadim=2,basin_only=False):
    '''
    lon,lat must be 2D mesh of coordinates
    '''
    dirlist = os.getcwd().split('/')[:-1]
    mydir=''
    for d in dirlist:
        mydir = mydir+d+'/'
#     print(mydir)
    if 'Abby.Hutson' in mydir:
        mydir='/home/Abby.Hutson/'
    # import datafile containing lon/lat points of basin polygons
    ds = xr.open_dataset(f'{mydir}GreatLakesandWa/greatlakes_basin_coords.nc')
    ont_points = ds.ont.data
    erie_points= ds.erie.data
    hur_points = ds.huron.data
    mich_points= ds.mich.data
    sup_points = ds.sup.data
    basin_points=ds.basin.data

    path      = mpth.Path(basin_points) # This is the entire basin, not just individual lakes
    basin_mask= path.contains_points(np.hstack((lon.flatten()[:,np.newaxis],lat.flatten()[:,np.newaxis])),
                                     radius=-0.1)
    rain_basin=basin_data(basin_mask,data,datadim)
    if basin_only is True:
        return rain_basin
    else:
        # Find locations that fall within the basin polygons
        path     = mpth.Path(mich_points) #<--create path using lat/lon data of polygon outline
        mich_mask= path.contains_points(np.hstack((lon.flatten()[:,np.newaxis],lat.flatten()[:,np.newaxis])),
                                        radius=-0.1)#<--find wrf points that fall within polygon
        path     = mpth.Path(ont_points)
        ont_mask = path.contains_points(np.hstack((lon.flatten()[:,np.newaxis],lat.flatten()[:,np.newaxis])),
                                        radius=-0.1)
        path     = mpth.Path(erie_points)
        erie_mask= path.contains_points(np.hstack((lon.flatten()[:,np.newaxis],lat.flatten()[:,np.newaxis])),
                                        radius=-0.1)
        path     = mpth.Path(hur_points)
        hur_mask = path.contains_points(np.hstack((lon.flatten()[:,np.newaxis],lat.flatten()[:,np.newaxis])),
                                        radius=-0.1)
        path     = mpth.Path(sup_points)
        sup_mask = path.contains_points(np.hstack((lon.flatten()[:,np.newaxis],lat.flatten()[:,np.newaxis])),
                                        radius=-0.1)

        # mask data within basins, set nan outside of the basins, append to year arrays
        rain_ont = basin_data(ont_mask,data)
        rain_erie= basin_data(erie_mask,data)
        rain_hur = basin_data(hur_mask,data)
        rain_mich= basin_data(mich_mask,data)
        rain_sup = basin_data(sup_mask,data)
        return rain_basin,rain_ont,rain_erie,rain_hur,rain_mich,rain_sup

def get_masked_prec_lake(ds, index, gpd_shapefile,varname='tp'):
    # splitting shapefile into individual polygons is more difficult than it needs to be

    # 0: erie, 1: huron, 2: michigan, 3: ontario, 4: superior (alphabetical)
    geometry_s = gpd_shapefile.geometry.apply(mapping)[index] # grab polygon for individual lake
    geometry   = json.dumps(geometry_s) # change single quotes to double quotes (dumb, right??)
    cropping_geoms = [geojson.loads(str(geometry))] # create geojson polygon object
    era_lake       =ds.rio.clip(geometries=cropping_geoms,
                                crs=gpd_shapefile.crs,
                                all_touched=True, drop=False) #clip based on geojson polygon

    # data is organized (month, hour, lat, lon)
#     tp_lake = np.sum(era_lake[varname].data,axis=1) THIS WAS ONLY SPECIFIC TO ERA DATA
    tp_lake = era_lake[varname].data
    return tp_lake

def mask_using_shapefile_polygon(dataset,lons,lats,varname='tp',lon_vname='lon',lat_vname='lat',whole_ds=False):
    '''
    This method can only be used when lons, lats are 1D arrays
    (i.e., when the grid is defined by regular lon/lat points
    as opposed to regular earth distance points)
    '''
    # get file directory correct based on file system
    dirlist = os.getcwd().split('/')[:-1]
    mydir=''
    for d in dirlist:
        mydir = mydir+d+'/'
    if 'Abby.Hutson' in mydir:
        mydir='/home/Abby.Hutson/'

    # print(np.nanmax(lons), np.nanmax(lats), np.nanmin(lons),np.nanmin(lats))
    lats_t,lons_t = np.meshgrid(lons,lats)
    lat, lon      = lats_t.T, lons_t.T

    shapefile = f'{mydir}GreatLakesandWa/greatlakes_subbasins/greatlakes_subbasins.shp'
    reader = Reader(shapefile)

    basins  = reader.records()
    lk_erie = [basin for basin in reader.records() if basin.attributes['merge']=='lk_erie'][0]
    lk_huron= [basin for basin in reader.records() if basin.attributes['merge']=='lk_huron'][0]
    lk_mich = [basin for basin in reader.records() if basin.attributes['merge']=='lk_mich'][0]
    lk_ont  = [basin for basin in reader.records() if basin.attributes['merge']=='lk_ont'][0]
    lk_sup  = [basin for basin in reader.records() if basin.attributes['merge']=='lk_sup'][0]

    # print(np.nanmax(dataset.t.data))
    with rasterio.Env(OGR_ENABLE_PARTIAL_REPROJECTION=True) as env:
        dataset.rio.set_spatial_dims(x_dim=lon_vname,y_dim=lat_vname,inplace=True)
        dataset.rio.write_crs('epsg:4326',inplace=True)
        basin   = gpd.read_file(shapefile, crs='epsg:4326')
        ds_basin=dataset.rio.clip(basin.geometry.apply(mapping), basin.crs,
                                      all_touched=True, drop=False)
    if whole_ds is True:
        # print(np.nanmax(ds_basin.t.data))
        return ds_basin
    else:
        rainbasin=ds_basin[varname].data

        rainerie= get_masked_prec_lake(dataset, 0,basin,varname=varname)
        rainhur = get_masked_prec_lake(dataset, 1,basin,varname=varname)
        rainmich= get_masked_prec_lake(dataset, 2,basin,varname=varname)
        rainont = get_masked_prec_lake(dataset, 3,basin,varname=varname)
        rainsup = get_masked_prec_lake(dataset, 4,basin,varname=varname)

        return rainbasin,rainont,rainerie,rainhur,rainmich,rainsup
