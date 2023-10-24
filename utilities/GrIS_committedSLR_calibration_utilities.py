#!python

import numpy as np
import xarray as xr

from scipy import stats
from scipy import integrate

import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import mascons

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from osgeo import ogr
import matplotlib.path as mplPath

import csv

import rasterio
from pyproj import Proj, transform, Transformer

import boto3

import warnings


def read_SERAC_obs(SERAC_obj_name):
    # Read dynamic dh from SERAC
    client = boto3.client('s3')
    obj = client.get_object(Bucket='dh-gapc', Key=SERAC_obj_name)
    data = obj['Body'].read().decode('utf-8').splitlines()

    records = csv.reader(data)
    headers = next(records)

    dh_dyn_obs = []

    #inProj  = Proj(init='epsg:32624')
    #outProj = Proj(init='epsg:3413')
    transformer = Transformer.from_crs('epsg:32624', 'epsg:3413')

    for row in records:
        dh_dyn_obs.append(dict.fromkeys(['PointID','x','y','dh0304', 'dh0405', 'dh0506', 'dh0607', \
                                                           'dh0708', 'dh0809', 'dh0910', 'dh1011', \
                                                           'dh1112', 'dh1213', 'dh1314', 'dh1415']))

        dh_dyn_obs[-1]['PointID'] = int(row[headers.index('PointID')])

        #xObs, yObs = transform(32624, 3413, float(row[headers.index('X')]), float(row[headers.index('Y')]))
        xObs, yObs = transformer.transform(float(row[headers.index('X')]), float(row[headers.index('Y')]))

        dh_dyn_obs[-1]['x'] = xObs
        dh_dyn_obs[-1]['y'] = yObs

        # Check FlgRangeH to decide whether to use DPSxxyy or DTSxxyy
        if float(row[headers.index('FlgRangeH')]) == 5:
            est_prefix = 'DPS' # ALPS ... more robust for large signals
            dh_sigma = 0.5 # m/yr
        elif float(row[headers.index('FlgRangeH')]) == -5:
            est_prefix = 'DTS' # polynomial ... more robust for small signals
            dh_sigma = 0.05 # m/yr
        else:
            est_prefix = ''
            dh_sigma = np.nan
        
        dh_dyn_obs[-1]['type'] = est_prefix
        
        SERAC_keys = [est_prefix + xxyy for xxyy in ['0304', '0405', '0506', '0607', '0708', '0809', \
                                                     '0910', '1011', '1112', '1213', '1314', '1415']]
        dh_dyn_obs_keys = [k.replace(est_prefix, 'dh') for k in SERAC_keys]
        
        for SERAC_key, dh_dyn_obs_keys in zip(SERAC_keys, dh_dyn_obs_keys):
            if row[headers.index(SERAC_key)]:
                dh_dyn_obs[-1][dh_dyn_obs_keys] = float(row[headers.index(SERAC_key)])
            else:
                dh_dyn_obs[-1][dh_dyn_obs_keys] = np.nan

        # Error is based on measurement type (DPS vs DTS)
        dh_dyn_obs[-1]['dh_sigma'] = dh_sigma
        
    print('total    number of obs: {:8d}'.format(len(dh_dyn_obs)))

    return dh_dyn_obs
    
def select_dh_obs(dhDYN_obs, startYear=2007, endYear=2015, \
                  allYears=True, method='sum', dhThresholdMax=None, dhThresholdMin=None, \
                  modelDomain=False, xModel=False, yModel=False, \
                  xMin=None, xMax=None, yMin=None, yMax=None): ##{{{
    
    startYear_string       = '{:4d}'.format(int(startYear  ))[2:]
    endYear_string         = '{:4d}'.format(int(endYear)    )[2:]
    
    # Sum annual dh over <start_year> - <end_year>
    # NOTE: Because SERAC dh obs cover the balance year (August-September), and the model output covers
    #       the calendar year (January-December), we throw in an extra year of dh from SERAC, which would
    #       capture dh over the Spring/Summer and give us a better match with the model output, which also
    #       captures Spring/Summer of the first year.

    dhDYN_obs_selected = []
    for iObs in range(len(dhDYN_obs)):
        dhOBS = list()
        dhOBS = [dhDYN_obs[iObs]['dh' + '{:4d}'.format(int(y-1))[2:] + '{:4d}'.format(int(y))[2:]] for y in range(int(startYear), int(endYear+1))]
        
        if method == 'sum':
            dhOBStotal = np.nansum(dhOBS)
        elif method == 'linear_fit':
            print('Interpolation not yet implemented')
            return None
    
        selectFlag = True
        if allYears:
            if np.any(np.isnan(dhOBS)):
                selectFlag = False
        if dhThresholdMax:
            if dhOBStotal > dhThresholdMax:
                selectFlag = False
        if dhThresholdMin:
            if dhOBStotal < dhThresholdMin:
                selectFlag = False
        if np.any(modelDomain):
            row = np.argmin(np.abs(yModel - dhDYN_obs[iObs]['y']))
            col = np.argmin(np.abs(xModel - dhDYN_obs[iObs]['x']))
            dh_mod_sample = modelDomain[row, col]
            if np.isnan(dh_mod_sample):
                selectFlag = False
        if xMin and xMax and yMin and yMax:
            if dhDYN_obs[iObs]['x'] <= xMin or dhDYN_obs[iObs]['x'] >= xMax or \
               dhDYN_obs[iObs]['y'] <= yMin or dhDYN_obs[iObs]['y'] >= yMax:
                selectFlag = False
                

        if selectFlag:
            dhDYN_obs_selected.append(dict.fromkeys(['PointID','x','y','dh'+startYear_string+endYear_string]))
            dhDYN_obs_selected[-1]['PointID'] = dhDYN_obs[iObs]['PointID']
            dhDYN_obs_selected[-1]['x'] = dhDYN_obs[iObs]['x']
            dhDYN_obs_selected[-1]['y'] = dhDYN_obs[iObs]['y']
            dhDYN_obs_selected[-1]['dh' + startYear_string + endYear_string] = dhOBStotal

            keys = ['dh' + '{:4d}'.format(int(y-1))[2:] + '{:4d}'.format(int(y))[2:] for y in range(int(startYear), int(endYear+1))]
            keys.append('dh_sigma')
            for key in keys:
                dhDYN_obs_selected[-1][key] = dhDYN_obs[iObs][key]
            
    return dhDYN_obs_selected
##}}}


def read_vel_error_tifs(year):
    # Read velocity
    vx_tif = 's3://dh-gapc/Velocity/nsidc-0478/greenland_vel_mosaic500_{:4.0f}_{:4.0f}_vx_v02.1.tif'.format(year, year+1)
    with rasterio.open(vx_tif) as ds:
        vx = ds.read(1)
        vx[vx == ds.nodata] = np.nan
        vv_transform = ds.transform
        vv_bounds = ds.bounds
    vy_tif = 's3://dh-gapc/Velocity/nsidc-0478/greenland_vel_mosaic500_{:4.0f}_{:4.0f}_vy_v02.1.tif'.format(year, year+1)
    with rasterio.open(vy_tif) as ds:
        vy = ds.read(1)
        vy[vy == ds.nodata] = np.nan
    vv = np.sqrt(vx**2 + vy**2)

    # Read error
    ex_tif = 's3://dh-gapc/Velocity/nsidc-0478/greenland_vel_mosaic500_{:4.0f}_{:4.0f}_ex_v02.1.tif'.format(year, year+1)
    with rasterio.open(ex_tif) as ds:
        ex = ds.read(1)
        ex[ex == ds.nodata] = np.nan
    ey_tif = 's3://dh-gapc/Velocity/nsidc-0478/greenland_vel_mosaic500_{:4.0f}_{:4.0f}_ey_v02.1.tif'.format(year, year+1)
    with rasterio.open(ey_tif) as ds:
        ey = ds.read(1)
        ey[ey == ds.nodata] = np.nan
    ev = np.sqrt(ex**2 + ey**2)
    
    # Write velocity
    
    # Write error
    
    return vx, vy, vv, ex, ey, ev, vv_transform, vv_bounds


# Model - observation dh differences
def grid_obs_dh(x_centers, y_centers, x_span, y_span, dh_dyn_obs_selected, startYear=None, endYear=None): ##{{{
    x_obs = np.array([d['x'] for d in dh_dyn_obs_selected])
    y_obs = np.array([d['y'] for d in dh_dyn_obs_selected])
    
    if startYear and endYear:
        startYear_string = '{:4d}'.format(int(startYear))[2:]
        endYear_string   = '{:4d}'.format(int(endYear)  )[2:]
        dh_obs = np.array([d['dh' + startYear_string + endYear_string] for d in dh_dyn_obs_selected])
        dh_sigma_obs = np.array([np.sqrt(endYear-startYear)*d['dh_sigma'] for d in dh_dyn_obs_selected])
    else:
        dh_obs = np.array([d['dh0715'] for d in dh_dyn_obs_selected])
        dh_sigma_obs = np.array([np.sqrt(8)*d['dh_sigma'] for d in dh_dyn_obs_selected])

    dh_obs_grid = np.nan * np.zeros( (len(y_centers), len(x_centers)) )
    dh_obs_sigma_grid = np.nan * np.zeros( (len(y_centers), len(x_centers)) )
    
    for col, x_center in enumerate(x_centers):
        for row, y_center in enumerate(y_centers):
            # obs
            idx1 = x_obs >= x_center-x_span/2
            idx2 = x_obs <= x_center+x_span/2
            idx3 = y_obs >= y_center-y_span/2 
            idx4 = y_obs <= y_center+y_span/2

            idx = np.logical_and( np.logical_and( np.logical_and( idx1, idx2), idx3), idx4)

            dh_obs_grid[row, col] = np.mean(dh_obs[idx])
            dh_obs_sigma_grid[row, col] = np.mean(dh_sigma_obs[idx])

    return dh_obs_grid, dh_obs_sigma_grid
##}}}

def grid_mod_dh(mod_name, x_centers, y_centers, x_span, y_span, var_name='dh_dynanom'): #{{{
    # Read model data
    gis_ds = xr.open_dataset('s3://' + mod_name, engine='zarr')

    x_mod = gis_ds['x'].data
    y_mod = gis_ds['y'].data

    z_mod = gis_ds[var_name]
    
    gis_ds.close()
    
    if var_name == 'dh_dynanom': z_mod = z_mod.fillna(np.nan).where(z_mod < 1000.)
    
    if len(z_mod.shape) == 3: z_mod = np.sum(z_mod.data, axis=0)
        
    z_mod_grid = np.nan * np.zeros( (len(y_centers), len(x_centers)) )
    
    for col, x_center in enumerate(x_centers):
        for row, y_center in enumerate(y_centers):
            
            idx1 = x_mod >= x_center-x_span/2
            idx2 = x_mod <= x_center+x_span/2
            idx3 = y_mod >= y_center-y_span/2 
            idx4 = y_mod <= y_center+y_span/2
            
            idx_col = np.where(np.logical_and(idx1, idx2))
            idx_row = np.where(np.logical_and(idx3, idx4))
            
            if len(idx_col[0]) > 0 and len(idx_row[0]) > 0:
                z_mod_grid[row, col] = np.nanmean(z_mod[np.min(idx_row):np.max(idx_row), \
                                                       np.min(idx_col):np.max(idx_col)])
            
    return z_mod_grid
##}}}


def grid_vel_obs(vel, vel_transform, x_centers, y_centers, x_span, y_span):
    height = vel.shape[0]
    width  = vel.shape[1]
    
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    xs, ys = rasterio.transform.xy(vel_transform, rows, cols)
    
    xs = np.array(xs); xs = xs[0,:]
    ys = np.array(ys); ys = ys[:,0]
    
    vel_obs_grid = np.nan * np.zeros( (len(y_centers), len(x_centers)) )

    for col, x_center in enumerate(x_centers):
        for row, y_center in enumerate(y_centers):
            
            idx1 = xs >= x_center-x_span/2
            idx2 = xs <= x_center+x_span/2
            idx3 = ys >= y_center-y_span/2 
            idx4 = ys <= y_center+y_span/2
            
            idx_col = np.where(np.logical_and(idx1, idx2))
            idx_row = np.where(np.logical_and(idx3, idx4))
            
            if len(idx_col[0]) > 0 and len(idx_row[0]) > 0:
               # Find number of valid observations
               nvalid = np.sum(~np.isnan(vel[np.min(idx_row):np.max(idx_row), np.min(idx_col):np.max(idx_col)]))
               if nvalid / ( (x_span/vel_transform[0]) * (y_span/-vel_transform[4]) ) >= 0.75:
                   vel_obs_grid[row, col] = np.nanmean(vel[np.min(idx_row):np.max(idx_row), \
                                                           np.min(idx_col):np.max(idx_col)])
               else:
                   vel_obs_grid[row, col] = np.nan
            
    return vel_obs_grid



## ----------------------------------------------------------------------------------------------------------------------------------------- ##
#def 
def calculate_residuals_and_grid(mod_name, x_centers, y_centers, grid_size, dh_obs_grid, dh_obs_sigma_grid, dh_mod_obs_sigma_multiplier):
    # Read model data
    gis_ds = xr.open_dataset('s3://' + mod_name, engine='zarr')

    x_mod = gis_ds['x'].data
    y_mod = gis_ds['y'].data

    dh_dynanom = gis_ds['dh_dynanom']
    
    gis_ds.close()
    
    dh_dynanom = dh_dynanom.fillna(np.nan).where(dh_dynanom < 1000.)
    
    dh_mod = np.sum(dh_dynanom, axis=0)
    dh_mod = dh_mod.fillna(np.nan).where(dh_mod < 1000.)
    
## ----------------------------------------------------------------------------------------------------------------------------------------- ##

def grid_models_and_calculate_residuals(mod_name, x_centers, y_centers, grid_size, obs_grid, obs_sigma_grid, mod_obs_sigma_multiplier, var_name='dh_dynanom'):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mod_grid = grid_mod_dh(mod_name, x_centers, y_centers, grid_size, grid_size, var_name=var_name)
    
    # Calculate residuals
    mod_obs_diff_grid = mod_grid - obs_grid
    
    # Flatten
    r = mod_obs_diff_grid.flatten()
    valid = ~np.isnan(r)
    r = r[valid]
    sigma = obs_sigma_grid.flatten()[valid]
    n = len(r)
        
    # Calculate score
    #s_j = np.exp( -(1./ 2.) * np.nansum( mod_obs_diff_grid**2 / (mod_obs_sigma_multiplier*obs_sigma_grid)**2 ) )
    s_j = np.exp( -(1./ 2.) * np.nansum( r**2 / (mod_obs_sigma_multiplier*sigma)**2 ) )
    #s_j = 1. / np.nansum( r**2 )

    # Calculate other statistics
    #rss = np.sqrt(np.sum(mod_obs_diff_grid[~np.isnan(mod_obs_diff_grid)]**2))
    
    return s_j, mod_obs_diff_grid

# Bayesian calibration
def gmsl_prior_posterior(gmsl, weights, xmin=0, xmax=15): ##{{{
    gkde_prior = stats.gaussian_kde(gmsl, bw_method='silverman')
    gkde_posterior = stats.gaussian_kde(gmsl, bw_method='silverman', weights=weights)

    x = np.linspace(xmin,xmax,10001)
    #x = np.linspace(0,np.max(gmsl),1001)
    kdepdf_prior = gkde_prior.evaluate(x)
    kdepdf_posterior = gkde_posterior.evaluate(x)

    return x, kdepdf_prior, kdepdf_posterior
##}}}

def stats_from_kde(x, y, thresholds=None):
    # Percentiles
    cdf = integrate.cumtrapz(y, x, initial=0)
    pcnt05 = x[np.where(cdf>=0.05)[0][0]]
    pcnt50 = x[np.where(cdf>=0.50)[0][0]]
    pcnt95 = x[np.where(cdf>=0.95)[0][0]]
    
    # MAP
    idx_map = np.argmax(y)
    MAP = x[idx_map]
    
    # Thresholds
    threshold_probabilities = None
    if thresholds is not None:
        threshold_probabilities = list()
        for threshold in thresholds:
            idx_threshold = np.where(x > threshold)[0][0]
            threshold_probabilities.append(100 * (1-cdf[idx_threshold]))

    return MAP, pcnt05, pcnt50, pcnt95, threshold_probabilities

def plot_GSFCmascons(lon_centers, lat_centers, cmwe_delta, min_lons, max_lons, min_lats, max_lats, ax=None, vmin=None, vmax=None, cmap='coolwarm_r'):
    polar_stereographic = ccrs.Stereographic(
        central_latitude=90.0,
        central_longitude=-45.0,
        false_easting=0.0,
        false_northing=0.0,
        true_scale_latitude=70.0,
        globe=ccrs.Globe('WGS84')
    )

    if vmin is None or vmax is None:
        min_mscns = np.min(cmwe_delta)
        max_mscns = np.max(cmwe_delta)
        vmax = np.max([np.abs(min_mscns), np.abs(max_mscns)])
        vmin = -vmax

    if ax is None:
        # Create figure and set projection:
        plt.figure(figsize=(5,5), dpi=300)
        ax = plt.axes(projection=polar_stereographic)
        
    ax.set_extent([-58, -27, 57, 84]) # Map bounds, [west, east, south, north]

    # In order to set a colormap for the mascon plot, we first plot a scatter of the mascons.
    # This will allow us to set colormap values for use with a second layer of "fill" data
    # representing the total extent of each mascon. These scatter points will be covered by
    # the "fill" images, and thus not visible in the final map.
    sc = ax.scatter(lon_centers, lat_centers, 1, c=cmwe_delta, zorder=0, transform=ccrs.PlateCarree(),
                         cmap=cmap, vmin=vmin, vmax=vmax)

    # Using the colormap info from above, draw each GIS mascon and fill with appropriate color
    cmap = plt.get_cmap(cmap)
    N_ints = 10
    for i in range(len(cmwe_delta)):
        x = np.append(np.linspace(min_lons[i], max_lons[i], N_ints),
                      np.linspace(max_lons[i], min_lons[i], N_ints))
        y = np.append(min_lats[i]*np.ones(N_ints), max_lats[i]*np.ones(N_ints))
        
        cmwe_delta_normalized = (cmwe_delta[i] - vmin) / (vmax - vmin)
        ax.fill(x, y, facecolor=cmap(cmwe_delta_normalized), edgecolor='none', zorder=5, transform=ccrs.PlateCarree())

    ax.coastlines(resolution='10m', zorder=7, linewidth=0.5)
    #ax.gridlines(zorder=8, linestyle=':', linewidth=0.5)
    
    sc.remove()
    
    return sc

def load_gscf_mascons():
    import numpy as np
    import xarray as xr
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import mascons

    gsfc = dict()
    
    file = open('../utilities/GSFC_highres_trends_GRACE_200701-201501_RL06v1.0_GSM-ICE6GD_GrIS.txt', 'r')
    lines = file.readlines()
    file.close()

    gsfc_cmwe_delta = list()
    gsfc_cmwe_delta_sigma = list()
    lat_centers = list()
    lon_centers = list()
    min_lons = list()
    max_lons = list()
    min_lats = list()
    max_lats = list()
    areas = list()
    labels = list()

    for line in lines:
        lat_centers.append(float(line.split(',')[0]))
        lon_centers.append(float(line.split(',')[1]))
        gsfc_cmwe_delta.append(float(line.split(',')[2]))
        gsfc_cmwe_delta_sigma.append(float(line.split(',')[3]))

        min_lats.append(float(line.split(',')[4]))
        max_lats.append(float(line.split(',')[5]))
        min_lons.append(float(line.split(',')[6]))
        max_lons.append(float(line.split(',')[7]))
        areas.append(float(line.split(',')[8]))
        labels.append(float(line.split(',')[9]))

    gsfc_cmwe_delta = np.array(gsfc_cmwe_delta)
    gsfc_cmwe_delta_sigma = np.array(gsfc_cmwe_delta_sigma)
    lat_centers = np.array(lat_centers)
    lon_centers = np.array(lon_centers)
    min_lats = np.array(min_lats)
    max_lats = np.array(max_lats)
    min_lons = np.array(min_lons)
    max_lons = np.array(max_lons)
    areas = np.array(areas)
    labels = np.array(labels)

    gsfc['cmwe_delta'] = np.array(gsfc_cmwe_delta)
    gsfc['cmwe_delta_sigma'] = np.array(gsfc_cmwe_delta_sigma)
    gsfc['lat_centers'] = np.array(lat_centers)
    gsfc['lon_centers'] = np.array(lon_centers)
    gsfc['min_lons'] = np.array(min_lons)
    gsfc['max_lons'] = np.array(max_lons)
    gsfc['min_lats'] = np.array(min_lats)
    gsfc['max_lats'] = np.array(max_lats)
    gsfc['areas'] = np.array(areas)
    gsfc['labels'] = np.array(labels)

    return gsfc
        
def plot_GSFCmascons_by_basin(lon_centers, lat_centers, gsfc_cmwe_delta, mascon_labels, ax=None, vmin=-5000, vmax=0, cmap='Reds_r', verbose=False):
    polar_stereographic = ccrs.Stereographic(
        central_latitude=90.0,
        central_longitude=-45.0,
        false_easting=0.0,
        false_northing=0.0,
        true_scale_latitude=70.0,
        globe=ccrs.Globe('WGS84')
    )
    

    xy_mascons = polar_stereographic.transform_points(ccrs.Geodetic(), lon_centers, lat_centers)

    gsfc_cmwe_delta_sums = dict()
    gsfc_cmwe_basin_rignot = dict()
    for basin_idx, basin_str in enumerate(['CW', 'NE', 'NO', 'NW', 'SE', 'SW']):
        basin_file = '../utilities/rignot_basins_' + basin_str + '.exp'

        x_basin = list()
        y_basin = list()
        xy_basin = list()
        with open(basin_file) as f:
            for i in range(5):
                next(f)
            for line in f:
                x_basin.append( float(line.split()[0]) )
                y_basin.append( float(line.split()[1]) )
                xy_basin.append( (float(line.split()[0]), float(line.split()[1]) ) )

        x_basin = np.array(x_basin)
        y_basin = np.array(y_basin)

        gsfc_cmwe_delta_sum = 0.
        for i in range(len(xy_mascons)):
            point = Point(xy_mascons[i][0], xy_mascons[i][1])
            polygon = Polygon(xy_basin)
            if polygon.contains(point):
                gsfc_cmwe_basin_rignot[mascon_labels[i]] = basin_str
                gsfc_cmwe_delta_sum += gsfc_cmwe_delta[i]
                
        gsfc_cmwe_delta_sums[basin_str] = gsfc_cmwe_delta_sum
        
        if verbose: print('{:2s} {: 4.1f}'.format(basin_str, gsfc_cmwe_delta_sum))
        
    sc = plot_values_by_basin(gsfc_cmwe_delta_sums, ax, cmap, vmin, vmax)

    return sc, gsfc_cmwe_basin_rignot


def plot_values_by_basin(gsfc_cmwe_delta_sums, ax, cmap, vmin, vmax):
    # Open the shapefile
    shapefile = ogr.Open('../utilities/IMBIE_Rignot_Basins/GRE_Basins_IMBIE2_v1.3.shp')

    # Get the layer
    layer = shapefile.GetLayer()

    # Iterate over the features in the layer
    x_basin = dict()
    y_basin = dict()
    x_means = dict()
    y_means = dict()
    for feature in layer:
        # Get the attributes of the feature
        attributes = feature.items()
        basin = attributes['SUBREGION1']
        if basin != 'ICE_CAP':
            # Get the geometry of the feature
            geometry = feature.GetGeometryRef()
            ring = geometry.GetGeometryRef(0)
            points = mplPath.Path(ring.GetPoints())
            x, y = zip(*points.to_polygons()[0])
            x_basin[basin] = x
            y_basin[basin] = y
            x_means[basin] = np.mean(x)
            y_means[basin] = np.mean(y)

    x_means_array = np.empty( len(gsfc_cmwe_delta_sums.keys()) )
    y_means_array = np.empty( len(gsfc_cmwe_delta_sums.keys()) )
    gsfc_cmwe_delta_sums_array = np.empty( len(gsfc_cmwe_delta_sums.keys()) )
    for i, basin in enumerate(gsfc_cmwe_delta_sums.keys()):
        x_means_array[i] = x_means[basin]
        y_means_array[i] = y_means[basin]
        gsfc_cmwe_delta_sums_array[i] = gsfc_cmwe_delta_sums[basin]

    sc = ax.scatter(x_means_array, y_means_array, 1, c=gsfc_cmwe_delta_sums_array, transform=ccrs.PlateCarree(),
                     cmap=cmap, vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(cmap)

    for i, basin in enumerate(gsfc_cmwe_delta_sums.keys()):
        gsfc_cmwe_delta_sums_array_normalized = (gsfc_cmwe_delta_sums_array[i] - vmin) / (vmax - vmin)
        ax.fill(x_basin[basin], y_basin[basin], facecolor=cmap(gsfc_cmwe_delta_sums_array_normalized), edgecolor='none', transform=ccrs.PlateCarree())

    ax.coastlines(resolution='10m', zorder=7, linewidth=0.5)
    #ax.gridlines(zorder=8, linestyle=':', linewidth=0.5)

    sc.remove()
    
    return sc