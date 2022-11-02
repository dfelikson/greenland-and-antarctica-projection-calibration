input_model_dir = './models';
output_netcdf_dir = './netcdfs';

vaf_cmmtt_csv_file = './proj-cmmtt-vaf-total-yearly.csv';
vaf_ctrl_csv_file  = './proj-ctrl-vaf-total-yearly.csv';

rho_ice = 917; % Hardcoded here but this is what is in md.materials.rho_ice for every ISSM simulation

times_hindcast = 2007:2020;
times_forecast = 2021:2101;
time = [times_hindcast times_forecast];

% read vaf csv
fid = fopen(vaf_cmmtt_csv_file, 'r');
vaf_cmmtt_header = textscan(strrep(fgetl(fid),'"',''),'%s','delimiter',',');
vaf_cmmtt_header = vaf_cmmtt_header{1};
vaf_cmmtt_cell = textscan(fid,'','delimiter',',');
fclose(fid);
fid = fopen(vaf_ctrl_csv_file, 'r');
vaf_ctrl_header = textscan(strrep(fgetl(fid),'"',''),'%s','delimiter',',');
vaf_ctrl_header = vaf_ctrl_header{1};
vaf_ctrl_cell = textscan(fid,'','delimiter',',');
fclose(fid);

d = dir(input_model_dir);
d_select = [];
for i = 1:length(d)
   if regexp(d(i).name, 'gris\.proj\.cmmtt\.A\d{4}')
      if isempty(d_select)
         d_select = d(i);
      else
         d_select(end+1) = d(i);
      end
   end
end

for ifile = 1:length(d_select) %%{{{

   spl = split(d_select(ifile).name, '.');
   ensembleID  = spl{4};
   
   fprintf(['processing ' d_select(ifile).name '\n']);

   output_thickness_netcdf_filename = [output_netcdf_dir '/gris.proj.cmmtt.' ensembleID '.thickness.nc'];
   output_velocity_netcdf_filename = [output_netcdf_dir '/gris.proj.cmmtt.' ensembleID '.velocity.nc'];
   output_MAF_netcdf_filename = [output_netcdf_dir '/gris.proj.cmmtt.' ensembleID '.MAF.nc'];
   if exist(output_thickness_netcdf_filename, 'file') & exist(output_velocity_netcdf_filename, 'file') & exist(output_MAF_netcdf_filename, 'file')
      fprintf(' -> all output netcdfs exist\n');
      continue
   end

   clear md

	% Check a few things to make sure we actually have to load the model, which can take a while
   % 1)  no thickness netcdf -> load md
   % 2) yes thickness netcdf
   %     no MAF netcdf +  no VAF in csv file -> load md
   % 3) yes thickness netcdf
   %     no MAF netcdf + yes VAF in csv file -> 
   input_model_filename = [input_model_dir '/' d_select(ifile).name];
   if ~exist(output_thickness_netcdf_filename, 'file')
      fprintf(' -> no thickness netcdf ... ');
      fprintf(['loading: ' input_model_filename '\n']);
      md = loadmodel(input_model_filename);
   end
   if ~exist(output_velocity_netcdf_filename, 'file')
      fprintf(' -> no velocity netcdf ... ');
      fprintf(['loading: ' input_model_filename '\n']);
      md = loadmodel(input_model_filename);
   end
   if exist(output_thickness_netcdf_filename, 'file') & ~exist(output_MAF_netcdf_filename, 'file') & ~any(strcmp(vaf_cmmtt_header,ensembleID))
      fprintf(' -> no MAF netcdf and no VAF in csv file ... ');
      fprintf(['loading: ' input_model_filename '\n']);
      md = loadmodel(input_model_filename);
   end
   if exist(output_thickness_netcdf_filename, 'file') & ~exist(output_MAF_netcdf_filename, 'file') & any(strcmp(vaf_cmmtt_header,ensembleID))
      fprintf(' -> no MAF netcdf and VAF found in csv file \n');
      maf = vaf_cmmtt_cell{find(strcmp(vaf_cmmtt_header,ensembleID))} * rho_ice;
   end

   if exist('md','var')
      times_all = [md.timestepping.start_time md.results.TransientSolution(:).time];
      idx_years = [];
      for year = [times_hindcast times_forecast]
         idx_years = [idx_years find(times_all == year)];
      end
      if times_all(1) ~= 2007
         disp('ERROR')
         return
      end
      thickness = [md.geometry.thickness md.results.TransientSolution(:).Thickness];
      velocity = [md.initialization.vel md.results.TransientSolution(:).Vel];
      mask = [md.mask.ice_levelset md.results.TransientSolution(:).MaskIceLevelset] < 0;
      thickness = thickness(:,idx_years);
      velocity = velocity(:,idx_years);
      mask = mask(:,idx_years);
      thickness(~mask) = NaN;
      velocity(~mask) = NaN;
      %vaf = [VolumeAboveFloatation(md) md.results.TransientSolution(:).IceVolumeAboveFloatation];
      vaf = [nan md.results.TransientSolution(:).IceVolumeAboveFloatation];
      maf = vaf(idx_years) * md.materials.rho_ice;

      md_mesh_x = md.mesh.x;
      md_mesh_y = md.mesh.y;
      md_mesh_elements = md.mesh.elements;
   end
      
   if ~exist(output_thickness_netcdf_filename, 'file')
      fprintf([' -> writing thickness netcdf: ' output_thickness_netcdf_filename '\n']);
      write_thickness_netcdf(md_mesh_elements, md_mesh_x, md_mesh_y, [times_hindcast times_forecast], thickness, output_thickness_netcdf_filename, 'lithk');
   end
   if ~exist(output_velocity_netcdf_filename, 'file')
      fprintf([' -> writing velocity netcdf: ' output_velocity_netcdf_filename '\n']);
		% NOTE: ISMIP6 did not define a velocity magnitude as part of standard output so, here, we use our own variable name (land_ice_surface_velocity) rather than an ISMIP6 variable name
      write_velocity_netcdf(md_mesh_elements, md_mesh_x, md_mesh_y, [times_hindcast times_forecast], velocity, output_velocity_netcdf_filename, 'land_ice_surface_velocity');
   end
   if ~exist(output_MAF_netcdf_filename, 'file')
      fprintf([' -> writing MAF netcdf: ' output_MAF_netcdf_filename '\n']);
      write_MAF_netcdf([times_hindcast times_forecast], maf, output_MAF_netcdf_filename, 'limnsw');
   end

end %%}}}

% ctrl
d = dir(input_model_dir);
d_select = [];
for i = 1:length(d)
   if regexp(d(i).name, 'gris\.proj\.ctrl\.A\d{4}')
      if isempty(d_select)
         d_select = d(i);
      else
         d_select(end+1) = d(i);
      end
   end
end

for ifile = 1:length(d_select) %%{{{

   spl = split(d_select(ifile).name, '.');
   ensembleID  = spl{4};
   
   fprintf(['processing ' d_select(ifile).name '\n']);

   output_thickness_netcdf_filename = [output_netcdf_dir '/gris.proj.ctrl.' ensembleID '.thickness.nc'];
   output_MAF_netcdf_filename = [output_netcdf_dir '/gris.proj.ctrl.' ensembleID '.MAF.nc'];
   if exist(output_thickness_netcdf_filename, 'file') & exist(output_velocity_netcdf_filename, 'file') & exist(output_MAF_netcdf_filename, 'file')
      fprintf(' -> all output netcdfs exist\n');
      continue
   end

   clear md

	% Check a few things to make sure we actually have to load the model, which can take a while
   % 1)  no thickness netcdf -> load md
   % 2) yes thickness netcdf
   %     no MAF netcdf +  no VAF in csv file -> load md
   % 3) yes thickness netcdf
   %     no MAF netcdf + yes VAF in csv file -> 
   input_model_filename = [input_model_dir '/' d_select(ifile).name];
   if ~exist(output_thickness_netcdf_filename, 'file')
      fprintf(' -> no thickness netcdf ... ');
      fprintf(['loading: ' input_model_filename '\n']);
      md = loadmodel(input_model_filename);
   end
   if exist(output_thickness_netcdf_filename, 'file') & ~exist(output_MAF_netcdf_filename, 'file') & ~any(strcmp(vaf_ctrl_header,ensembleID))
      fprintf(' -> no MAF netcdf and no VAF in csv file ... ');
      fprintf(['loading: ' input_model_filename '\n']);
      md = loadmodel(input_model_filename);
   end
   if exist(output_thickness_netcdf_filename, 'file') & ~exist(output_MAF_netcdf_filename, 'file') & any(strcmp(vaf_ctrl_header,ensembleID))
      fprintf(' -> no MAF netcdf and VAF found in csv file \n');
      vaf = vaf_ctrl_cell{find(strcmp(vaf_ctrl_header,ensembleID))} * rho_ice;
   end

   if exist('md','var')
      times_all = [md.timestepping.start_time md.results.TransientSolution(:).time];
      idx_years = [];
      for year = [times_hindcast times_forecast];
         idx_years = [idx_years find(times_all == year)];
      end
      if times_all(1) ~= 2007
         disp('ERROR')
         return
      end
      thickness = [md.geometry.thickness md.results.TransientSolution(:).Thickness];
      mask = repmat(md.mask.ice_levelset, 1, size(thickness,2)) < 0;
      thickness = thickness(:,idx_years);
      mask = mask(:,idx_years);
      thickness(~mask) = NaN;
      %vaf = [VolumeAboveFloatation(md) md.results.TransientSolution(:).IceVolumeAboveFloatation];
      vaf = [nan md.results.TransientSolution(:).IceVolumeAboveFloatation];
      maf = vaf(idx_years) * md.materials.rho_ice;

      md_mesh_x = md.mesh.x;
      md_mesh_y = md.mesh.y;
      md_mesh_elements = md.mesh.elements;
   end
      
   if ~exist(output_thickness_netcdf_filename, 'file')
      fprintf([' -> writing thickness netcdf: ' output_thickness_netcdf_filename '\n']);
      write_thickness_netcdf(md_mesh_elements, md_mesh_x, md_mesh_y, [times_hindcast times_forecast], thickness, output_thickness_netcdf_filename, 'lithk');
   end
   if ~exist(output_MAF_netcdf_filename, 'file')
      fprintf([' -> writing MAF netcdf: ' output_MAF_netcdf_filename '\n']);
      write_MAF_netcdf([times_hindcast times_forecast], vaf, output_MAF_netcdf_filename, 'limnsw');
   end

end %%}}}

