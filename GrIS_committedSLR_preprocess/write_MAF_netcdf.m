function write_MAF_netcdf(time, MAF, output_netcdf_filename, variable_name)

   % time
   nccreate(output_netcdf_filename,'time','Dimensions',{'time',length(time)})
   ncwrite(output_netcdf_filename,'time',time)
   ncwriteatt(output_netcdf_filename,'time','standard_name','time')
   ncwriteatt(output_netcdf_filename,'time','long_name','time (year)')

   % MAF
   nccreate(output_netcdf_filename,variable_name,'FillValue',9.96921e36,'Dimensions',{'time',length(time)})
   ncwrite(output_netcdf_filename,variable_name, MAF)
   ncwriteatt(output_netcdf_filename,variable_name,'grid_mapping','Polar_Stereographic')
   ncwriteatt(output_netcdf_filename,variable_name,'standard_name','land_ice_mass_not_displacing_sea_water')
   ncwriteatt(output_netcdf_filename,variable_name,'units','kg')

end % main function

