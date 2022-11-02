function write_velocity_netcdf(md_mesh_elements, md_mesh_x, md_mesh_y, ...
      time, velocity, output_netcdf_filename, variable_name)

   % for interpolation
   xgrid = -720000:1000:960000; % projected x coordinates spaced 1 km apart
   ygrid = -3450000:1000:-570000; % projected y coordinates spaced 1 km apart
   
   % interpolate onto grid
   velocity_grid = zeros(length(ygrid), length(xgrid), length(time));
   for j = 1:length(time)
       [~, velocity_grid(:,:,j)] = evalc('InterpFromMeshToGrid(md_mesh_elements, md_mesh_x, md_mesh_y, velocity(:,j), xgrid, ygrid, NaN);');
   end

   % % interpolate dh_smb onto grid
   % if size(dh_smb,2) > 1
   %    dh_smb_grid = zeros(length(ygrid), length(xgrid), length(times_hindcast));
   %    for j = 1:length(times_hindcast)
   %        [~, dh_smb_grid(:,:,j)] = evalc('InterpFromMeshToGrid(md_mesh_elements, md_mesh_x, md_mesh_y, dh_smb(:,j), xgrid, ygrid, NaN);');
   %    end
   % else
   %    dh_smb_grid = zeros(length(ygrid), length(xgrid), 1);
   %    [~, dh_smb_grid(:,:)] = evalc('InterpFromMeshToGrid(md_mesh_elements, md_mesh_x, md_mesh_y, dh_smb, xgrid, ygrid, NaN);');
   % end

   % coordinates and times {{{
   % x
   nccreate(output_netcdf_filename,'x','Dimensions',{'x',length(xgrid)})
   ncwrite(output_netcdf_filename,'x',xgrid)
   ncwriteatt(output_netcdf_filename,'x','units','m')
   ncwriteatt(output_netcdf_filename,'x','standard_name','projection_x_coordinate')
   ncwriteatt(output_netcdf_filename,'x','long_name','x coordinate of projection')

   % y
   nccreate(output_netcdf_filename,'y','Dimensions',{'y',length(ygrid)})
   ncwrite(output_netcdf_filename,'y',ygrid)
   ncwriteatt(output_netcdf_filename,'y','units','m')
   ncwriteatt(output_netcdf_filename,'y','standard_name','projection_y_coordinate')
   ncwriteatt(output_netcdf_filename,'y','long_name','y coordinate of projection')
   
   % time
   nccreate(output_netcdf_filename,'time','Dimensions',{'time',length(time)})
   ncwrite(output_netcdf_filename,'time',time)
   ncwriteatt(output_netcdf_filename,'time','standard_name','time')
   ncwriteatt(output_netcdf_filename,'time','long_name','time (year)')
   %}}}

   % projection {{{
   nccreate(output_netcdf_filename,'Polar_Stereographic','Datatype','char')
   ncwriteatt(output_netcdf_filename,'Polar_Stereographic','ellipsoid','WGS84')
   ncwriteatt(output_netcdf_filename,'Polar_Stereographic','false_easting',0.0)
   ncwriteatt(output_netcdf_filename,'Polar_Stereographic','false_northing',0.0)
   ncwriteatt(output_netcdf_filename,'Polar_Stereographic','grid_mapping_name','polar_stereographic')
   ncwriteatt(output_netcdf_filename,'Polar_Stereographic','latitude_of_projection_origin',90.0)
   ncwriteatt(output_netcdf_filename,'Polar_Stereographic','standard_parallel',70.0)
   ncwriteatt(output_netcdf_filename,'Polar_Stereographic','straight_vertical_longitude_from_pole',-45.0)
   %}}}

   % velocity
   nccreate(output_netcdf_filename,variable_name,'FillValue',9.96921e36,'Dimensions',{'y',length(ygrid),'x',length(xgrid),'time',length(time)})
   ncwrite(output_netcdf_filename,variable_name, velocity_grid)
   ncwriteatt(output_netcdf_filename,variable_name,'grid_mapping','Polar_Stereographic')
   ncwriteatt(output_netcdf_filename,variable_name,'standard_name','land_ice_surface_velocity')
   ncwriteatt(output_netcdf_filename,variable_name,'units','m/s')

   % % dh_smb
   % if length(size(dh_smb_cmmtt_grid)) > 2
   %    nccreate(output_netcdf_filename,'dh_smb_cmmtt','FillValue',9.96921e36,'Dimensions',{'x',length(xgrid),'y',length(ygrid),'time_hindcast',length(times_hindcast)})
   %    ncwrite(output_netcdf_filename,'dh_smb_cmmtt',permute(dh_smb_cmmtt_grid,[2 1 3]))
   % else
   %    nccreate(output_netcdf_filename,'dh_smb_cmmtt','FillValue',9.96921e36,'Dimensions',{'x',length(xgrid),'y',length(ygrid),'time_1year',1})
   %    ncwrite(output_netcdf_filename,'dh_smb_cmmtt',permute(dh_smb_cmmtt_grid,[2 1]))
   % end
   % ncwriteatt(output_netcdf_filename,'dh_smb_cmmtt','grid_Polar_Stereographic','mapping')
   % ncwriteatt(output_netcdf_filename,'dh_smb_cmmtt','standard_name','dh_smb_cmmtt')
   % ncwriteatt(output_netcdf_filename,'dh_smb_cmmtt','long_name','committed smb dh for the previous year')
   % ncwriteatt(output_netcdf_filename,'dh_smb_cmmtt','units','m')

   % % dVAF
   % nccreate(output_netcdf_filename,'dVAF_cmmtt','FillValue',9.96921e36,'Dimensions',{'time',length(time)})
   % ncwrite(output_netcdf_filename,'dVAF_cmmtt',dVAF_cmmtt)
   % ncwriteatt(output_netcdf_filename,'dVAF_cmmtt','grid_Polar_Stereographic','mapping')
   % ncwriteatt(output_netcdf_filename,'dVAF_cmmtt','standard_name','dVAF_cmmtt')
   % ncwriteatt(output_netcdf_filename,'dVAF_cmmtt','long_name','committed dVAF for the previous year')
   % ncwriteatt(output_netcdf_filename,'dVAF_cmmtt','units','m3')

end % main function

