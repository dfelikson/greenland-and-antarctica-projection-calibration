#!/usr/bin/env python

# Imports ##{{{
import os, sys, glob, time

import progressbar

from netCDF4 import Dataset

import numpy as np
import scipy
from scipy import stats

from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

sys.path.append('/Users/dfelikso/Research/Software/ScriptsAndUtilities/pythonModules')
import raster
import shapefile_utils
import bilinear_interpolate

# GRACE mascon comparison
sys.path.append('/Users/dfelikso/Research/Software/ScriptsAndUtilities/pythonModules/Grace_GIS_Mascon')
import mascons

import GrIS_committedSLR_calibration_utilities
##}}}

# Setup ##{{{
model_netcdf_dir = './netcdfs'

start_year = 2007.
end_year = 2015.
start_year2 = start_year - 2000.
end_year2 = end_year - 2000.

orange = [c/255 for c in [255, 150, 0]]
darkorange = [c/255 for c in [255, 125, 0]]

# plot axis
gris_axis_km = [-600,  800, -3350,  -750]
jak_axis_km  = [-195, -150, -2290, -2261]
hel_axis_km  = [ 265,  320, -2590, -2520]
nw_glacier_axis_km = [-408, -344, -1489, -1432]

# Read model files
cmmtt_files = glob.glob(model_netcdf_dir + '/gris.proj.cmmtt.A????.thickness.nc')
velocity_files = glob.glob(model_netcdf_dir + '/gris.proj.cmmtt.A????.velocity.nc')

# coordinates
#cmmtt_file = cmmtt_files[0]
#cmmtt_ds = Dataset(cmmtt_file)

x = np.arange(-720000,961000,1000) #cmmtt_ds.variables['x'][:]
y = np.arange(-3450000,-571000,1000) #cmmtt_ds.variables['y'][:]
xstep  = x[ 1] - x[0]
ystep  = y[ 1] - y[0]
left   = x[ 0] - xstep/2
right  = x[-1] + xstep/2
bottom = y[ 0] - ystep/2
top    = y[-1] + ystep/2
extent = [m/1000 for m in [left, right, bottom, top]]

# time
start_idx = 2007 #np.where(cmmtt_ds['time'][:] == start_year)[0][0]
end_idx = 2014 #np.where(cmmtt_ds['time'][:] == end_year)[0][0]

#cmmtt_ds.close()

# smb baseline
#smb6089_ds = Dataset(model_netcdf_dir + '/gris.smb6089.nc')
#smb6089 = smb6089_ds['smb'][:,:].transpose()
#smb6089_ds.close()
##}}}

# Read model files
overwrite = False

dh_anomaly_sum = None
dh_dynAnom_sum = None

for iFile, cmmtt_file in enumerate(cmmtt_files): ##{{{
   ensembleID = os.path.basename(cmmtt_file).split('.')[3]
   # netcdfFile_dhanomaly = model_netcdf_dir + \
   #    '/gris.proj.{:4.0f}-{:4.0f}.{:s}.dhanomaly.nc'.format(start_year, end_year, ensembleID)
   # netcdfFile_dhdyn = model_netcdf_dir + \
   #    '/gris.proj.{:4.0f}-{:4.0f}.{:s}.dhdynamic.nc'.format(start_year, end_year, ensembleID)
   netcdfFile_dhdynAnom = model_netcdf_dir + \
      '/gris.proj.{:4.0f}-{:4.0f}.{:s}.dhdynAnom6079.nc'.format(start_year, end_year, ensembleID)         
      
   #if os.path.exists(netcdfFile_dhanomaly) and os.path.exists(netcdfFile_dhdyn) and not overwrite:
   if os.path.exists(netcdfFile_dhdynAnom) and not overwrite:
      print('skipping ensemble member {:s}: {:s}'.format(ensembleID, os.path.basename(cmmtt_file)))
   else:
      print('loading ensemble member {:s}: {:s}'.format(ensembleID, os.path.basename(cmmtt_file)))
      
      # cmmtt
      cmmtt_ds = Dataset(cmmtt_file)
      start_idx = np.where(cmmtt_ds['time'][:] == start_year)[0][0]
      end_idx = np.where(cmmtt_ds['time'][:] == end_year)[0][0]
      lithk = cmmtt_ds['lithk'][start_idx:end_idx+1,:,:]
      dh_cmmtt = np.diff(lithk, axis=0).filled(np.nan).transpose(0,2,1)

      # # ctrl
      # ctrl_file = cmmtt_file.replace('cmmtt','ctrl')
      # if not os.path.exists(ctrl_file):
      #    print(' -> WARNING: File {:s} does not exist!'.format(os.path.basename(ctrl_file)))
      # else:
      #    print('                         {:s}'.format(os.path.basename(ctrl_file)))
      #    ctrl_ds = Dataset(ctrl_file)
      #    lithk_ctrl = ctrl_ds['lithk'][start_idx:end_idx+1,:,:]
      #    dh_ctrl  = np.diff(lithk_ctrl, axis=0).filled(np.nan).transpose(0,2,1)

      # smb_anom
      smb_anom_file = model_netcdf_dir + '/gris.smb_anom6079.' + ensembleID + '.nc'
      if not os.path.exists(smb_anom_file):
         print(' -> WARNING: File {:s} does not exist!'.format(os.path.basename(smb_anom_file)))
      else:
         print('                               {:s}'.format(os.path.basename(smb_anom_file)))
         smb_anom_ds = Dataset(smb_anom_file)
         smb_anom = smb_anom_ds['smb_anom'][:,:].transpose()
         #smb_anom[smb_anom==-9999.] = np.nan
         dh_smb_anom = np.repeat(smb_anom[np.newaxis, :, :], end_year-start_year, axis=0)

      # # dh_anomaly
      # dh_anomaly = dh_cmmtt - dh_ctrl
      # print(' -> writing dhanomaly ncfile: ' + os.path.basename(netcdfFile_dhanomaly))
      # if os.path.exists(netcdfFile_dhanomaly):
      #    os.remove(netcdfFile_dhanomaly)
      # ncfile = Dataset(netcdfFile_dhanomaly,'w')
      # ncfile.createDimension('t', dh_ctrl.shape[0])
      # ncfile.createDimension('x', dh_ctrl.shape[2])
      # ncfile.createDimension('y', dh_ctrl.shape[1])
      # # t/x/y coordinates
      # t_var = ncfile.createVariable('t','f4',('t',))
      # t_var[:] = np.arange(int(start_year), int(end_year))
      # x_var = ncfile.createVariable('x','f4',('x',))
      # x_var[:] = cmmtt_ds.variables['x'][:]
      # y_var = ncfile.createVariable('y','f4',('y',))
      # y_var[:] = cmmtt_ds.variables['y'][:]
      # # dh_anomaly
      # dh_anomaly_var = ncfile.createVariable('dh_anomaly','f4',('t', 'y', 'x'))
      # dh_anomaly_var[:,:,:] = dh_anomaly#.transpose(0,2,1)
      # ncfile.close()
      # 
      # if dh_anomaly_sum is None:
      #    dh_anomaly_sum = dh_anomaly
      # else:
      #    dh_anomaly_sum = dh_anomaly_sum + dh_anomaly

      # # dh_dynamic
      # dh_dynamic = dh_cmmtt - dh_smb
      # print(' -> writing dhdynamic ncfile: ' + os.path.basename(netcdfFile_dhdyn))
      # if os.path.exists(netcdfFile_dhdyn):
      #    os.remove(netcdfFile_dhdyn)
      # ncfile = Dataset(netcdfFile_dhdyn,'w')
      # ncfile.createDimension('t', dh_ctrl.shape[0])
      # ncfile.createDimension('x', dh_ctrl.shape[2])
      # ncfile.createDimension('y', dh_ctrl.shape[1])
      # # t/x/y coordinates
      # t_var = ncfile.createVariable('t','f4',('t',))
      # t_var[:] = np.arange(int(start_year), int(end_year))
      # x_var = ncfile.createVariable('x','f4',('x',))
      # x_var[:] = cmmtt_ds.variables['x'][:]
      # y_var = ncfile.createVariable('y','f4',('y',))
      # y_var[:] = cmmtt_ds.variables['y'][:]
      # # dh_dynamic
      # dh_dynamic_var = ncfile.createVariable('dh_dynamic','f4',('t', 'y', 'x'))
      # dh_dynamic_var[:,:,:] = dh_dynamic#.transpose(0,2,1)
      # ncfile.close()

      # dh_dynamic_anomaly
      print('')
      print(cmmtt_file)
      print(smb_anom_file)
      dh_dynAnom = dh_cmmtt - dh_smb_anom

      # ## DEBUG
      # fig, ax = plt.subplots(1,3,figsize=(15,8), sharex=True, sharey=True)

      # idx = 0
      # vmin = -20
      # vmax = +20
      # im_dh = ax[0].imshow(np.sum(dh_cmmtt, axis=0), \
      #            cmap='RdBu', vmin=vmin, vmax=vmax, origin='lower', \
      #            aspect='equal')

      # im_dh = ax[1].imshow(np.sum(dh_smb_anom, axis=0), \
      #            cmap='RdBu', vmin=vmin, vmax=vmax, origin='lower', \
      #            aspect='equal')

      # im_dh = ax[2].imshow(np.sum(dh_dynAnom, axis=0), \
      #            cmap='RdBu', vmin=vmin, vmax=vmax, origin='lower', \
      #            aspect='equal')

      # plt.show()
      # import pdb; pdb.set_trace()
      # ## DEBUG

      print(' -> writing dhdynAnom ncfile:  ' + os.path.basename(netcdfFile_dhdynAnom))
      if os.path.exists(netcdfFile_dhdynAnom):
         os.remove(netcdfFile_dhdynAnom)
      ncfile = Dataset(netcdfFile_dhdynAnom,'w')
      ncfile.createDimension('t', dh_cmmtt.shape[0])
      ncfile.createDimension('x', dh_cmmtt.shape[2])
      ncfile.createDimension('y', dh_cmmtt.shape[1])
      # t/x/y coordinates
      t_var = ncfile.createVariable('t','f4',('t',))
      t_var[:] = np.arange(int(start_year), int(end_year))
      x_var = ncfile.createVariable('x','f4',('x',))
      x_var[:] = cmmtt_ds.variables['x'][:]
      y_var = ncfile.createVariable('y','f4',('y',))
      y_var[:] = cmmtt_ds.variables['y'][:]
      # dh_dynamic
      dh_dynamic_var = ncfile.createVariable('dh_dynAnom','f4',('t', 'y', 'x'))
      dh_dynamic_var[:,:,:] = dh_dynAnom#.transpose(0,2,1)
      ncfile.close()
      
      if dh_dynAnom_sum is None:
         dh_dynAnom_sum = dh_dynAnom
      else:
         dh_dynAnom_sum = dh_dynAnom_sum + dh_dynAnom
      
      cmmtt_ds.close()
      # ctrl_ds.close()
      # smb_ds.close()
      print('')
##}}}

print('')
for iFile, velocity_file in enumerate(velocity_files): ##{{{
   ensembleID = os.path.basename(velocity_file).split('.')[3]
   netcdfFile_dv = model_netcdf_dir + \
      '/gris.proj.{:4.0f}-{:4.0f}.{:s}.dv.nc'.format(start_year, end_year, ensembleID)         
      
   #if os.path.exists(netcdfFile_dhanomaly) and os.path.exists(netcdfFile_dhdyn) and not overwrite:
   if os.path.exists(netcdfFile_dv) and not overwrite:
      print('skipping ensemble member {:s}: {:s}'.format(ensembleID, os.path.basename(velocity_file)))
   else:
      print('loading ensemble member {:s}: {:s}'.format(ensembleID, os.path.basename(velocity_file)))
      
      # cmmtt
      v_ds = Dataset(velocity_file)
      start_idx = np.where(v_ds['time'][:] == start_year)[0][0]
      end_idx = np.where(v_ds['time'][:] == end_year)[0][0]
      velsurf = v_ds['land_ice_surface_velocity'][start_idx:end_idx+1,:,:]
      dv = np.diff(velsurf, axis=0).filled(np.nan).transpose(0,2,1)

      print(' -> writing dv ncfile:  ' + os.path.basename(netcdfFile_dv))
      if os.path.exists(netcdfFile_dv):
         os.remove(netcdfFile_dv)
      ncfile = Dataset(netcdfFile_dv,'w')
      ncfile.createDimension('t', dv.shape[0])
      ncfile.createDimension('x', dv.shape[2])
      ncfile.createDimension('y', dv.shape[1])
      # t/x/y coordinates
      t_var = ncfile.createVariable('t','f4',('t',))
      t_var[:] = np.arange(int(start_year), int(end_year))
      x_var = ncfile.createVariable('x','f4',('x',))
      x_var[:] = v_ds.variables['x'][:]
      y_var = ncfile.createVariable('y','f4',('y',))
      y_var[:] = v_ds.variables['y'][:]
      # dv
      dv_var = ncfile.createVariable('dv','f4',('t', 'y', 'x'))
      dv_var[:,:,:] = dv#.transpose(0,2,1)
      ncfile.close()
      
      v_ds.close()
      print('')
##}}}

# Calculate ensemble means
# dh_anomaly_mean_filename = model_netcdf_dir + \
#       '/gris.proj.{:4.0f}-{:4.0f}.dhanomalyMean.nc'.format(start_year, end_year, ensembleID)
# if os.path.exists(dh_anomaly_mean_filename):
#    print(' -> dh_anomaly_mean ncfile already exists')
#    ncfile = Dataset(dh_anomaly_mean_filename, 'r')
#    dh_anomaly_mean = ncfile['dh_anomalyMean'][:,:]
#    ncfile.close()
# else:
#    print(' -> writing dh_anomaly_mean ncfile')
#    dh_anomaly_mean = dh_anomaly_sum / len(cmmtt_files)
#    ncfile = Dataset(dh_anomaly_mean_filename, 'w')
#    ncfile.createDimension('t', dh_ctrl.shape[0])
#    ncfile.createDimension('x', dh_ctrl.shape[2])
#    ncfile.createDimension('y', dh_ctrl.shape[1])
#    # t/x/y coordinates
#    t_var = ncfile.createVariable('t','f4',('t',))
#    t_var[:] = np.arange(int(start_year), int(end_year))
#    x_var = ncfile.createVariable('x','f4',('x',))
#    x_var[:] = x
#    y_var = ncfile.createVariable('y','f4',('y',))
#    y_var[:] = y
#    # dh_dynamic
#    dh_dynamic_var = ncfile.createVariable('dh_anomalyMean','f4',('t', 'y', 'x'))
#    dh_dynamic_var[:,:,:] = dh_anomaly_mean#.transpose(0,2,1)
#    ncfile.close()

# dh_dynAnom_mean_filename = model_netcdf_dir + \
#       '/gris.proj.{:4.0f}-{:4.0f}.dhdynAnomMean.nc'.format(start_year, end_year, ensembleID)
# if os.path.exists(dh_dynAnom_mean_filename):
#    print(' -> dh_dynAnom_mean ncfile already exists')
#    ncfile = Dataset(dh_dynAnom_mean_filename, 'r')
#    dh_dynAnom_mean = ncfile['dh_dynAnomMean'][:,:]
#    ncfile.close()
# else:
#    print(' -> writing dh_dynAnom_mean ncfile')
#    dh_dynAnom_mean = dh_dynAnom_sum / len(cmmtt_files)
#    ncfile = Dataset(dh_dynAnom_mean_filename, 'w')
#    ncfile.createDimension('t', dh_ctrl.shape[0])
#    ncfile.createDimension('x', dh_ctrl.shape[2])
#    ncfile.createDimension('y', dh_ctrl.shape[1])
#    # t/x/y coordinates
#    t_var = ncfile.createVariable('t','f4',('t',))
#    t_var[:] = np.arange(int(start_year), int(end_year))
#    x_var = ncfile.createVariable('x','f4',('x',))
#    x_var[:] = x
#    y_var = ncfile.createVariable('y','f4',('y',))
#    y_var[:] = y
#    # dh_dynamic
#    dh_dynamic_var = ncfile.createVariable('dh_dynAnomMean','f4',('t', 'y', 'x'))
#    dh_dynamic_var[:,:,:] = dh_dynAnom_mean#.transpose(0,2,1)
#    ncfile.close()
# 
