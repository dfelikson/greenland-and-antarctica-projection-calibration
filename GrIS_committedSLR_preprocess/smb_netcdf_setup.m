% read in LHS table and initialise the model ready to run
% committed sea level projections
% by altering friction, viscosity and smb as indicated
% by the LHS factors

% read in ensemble table ("ensemble-128-4p.tab")
% This table provides the LHS values (between 0 and 1) rather than the
% factors used to vary the parameters, because those will be changed in here.
%lhs = tdfread('/Users/inias/Documents/science/committedSLR/ensemble-generation/ensemble-128-4p-noID.tab',' ');
file = fopen('ensemble-128-4p-noID.tab','r');
lhs_cell = textscan(file, '%f %f %f %f', 'EndOfLine', '\r\n', 'headerlines', 1);
fclose(file);
lhs.fric = lhs_cell{1};
lhs.visc = lhs_cell{2};
lhs.smb1 = lhs_cell{3};
lhs.smb2 = lhs_cell{4};

%file = fopen('/Users/inias/Documents/science/committedSLR/ensemble-generation/IDs.txt','r');
file = fopen('IDs.txt','r');
IDs = textscan(file,'%s', 'EndOfLine', '\r\n');
lhs.IDs = IDs{:, 1};
fclose(file);


% Load model
md = loadmodel('./models/gris.proj.cmmtt.A1091');

% SMB: load in 2001-2015 period -- perturbations will be added to this smb
% load in RACMO SMB
ncsmb = './data/smb_rec.1958-2017.BN_RACMO2.3p2_FGRN055_GrIS.MM.nc';
ncrunoff = './data/runoff.1958-2017.BN_RACMO2.3p2_FGRN055_GrIS.MM.nc';
x    = ncread(ncrunoff,'lon'); % despite name in netcdf, it isn't actually lat and lon, but projected x y
y    = ncread(ncrunoff,'lat');
smb_0115   = (sum(ncread(ncsmb,'SMB_rec',[1 1 517],[Inf Inf 180]),3)/15)'; % sum monthly data and divide by number of years
% Interpolate onto model mesh
smb_0115=InterpFromGridToMesh(x,y,smb_0115,md.mesh.x,md.mesh.y,0);
% convert mm/yr to m/yr of ice equivalent
smb_0115=(smb_0115*md.materials.rho_freshwater/md.materials.rho_ice)/1000; %to get m/yr ice equivalent

% SMB: find seasonality (for 1960-1989 reference period) %%{{{
% NOTE: This is used to calculate the SMB that's applied as forcing in the model runs.
nyears = 30;
nmonths = nyears*12;
% monthly data from 30 years from Jan 1960 - Dec 1989:
smbmonth   = ncread(ncsmb,'SMB_rec',[1 1 25],[Inf Inf nmonths]); % monthly smb in mm WE / month
% sequence of each Januaries
jans = 1:12:(nyears*12);
% mean monthly data for 1960-89
% converted from monthly mm WE totals to m/yr ice equivalent
% approximate that all months are 1/12 of a year
smb.jan = ((sum(smbmonth(:,:,jans   ),3)*12/nyears)*md.materials.rho_freshwater/md.materials.rho_ice)'/1000;
smb.feb = ((sum(smbmonth(:,:,jans+ 1),3)*12/nyears)*md.materials.rho_freshwater/md.materials.rho_ice)'/1000;
smb.mar = ((sum(smbmonth(:,:,jans+ 2),3)*12/nyears)*md.materials.rho_freshwater/md.materials.rho_ice)'/1000;
smb.apr = ((sum(smbmonth(:,:,jans+ 3),3)*12/nyears)*md.materials.rho_freshwater/md.materials.rho_ice)'/1000;
smb.may = ((sum(smbmonth(:,:,jans+ 4),3)*12/nyears)*md.materials.rho_freshwater/md.materials.rho_ice)'/1000;
smb.jun = ((sum(smbmonth(:,:,jans+ 5),3)*12/nyears)*md.materials.rho_freshwater/md.materials.rho_ice)'/1000;
smb.jul = ((sum(smbmonth(:,:,jans+ 6),3)*12/nyears)*md.materials.rho_freshwater/md.materials.rho_ice)'/1000;
smb.aug = ((sum(smbmonth(:,:,jans+ 7),3)*12/nyears)*md.materials.rho_freshwater/md.materials.rho_ice)'/1000;
smb.sep = ((sum(smbmonth(:,:,jans+ 8),3)*12/nyears)*md.materials.rho_freshwater/md.materials.rho_ice)'/1000;
smb.oct = ((sum(smbmonth(:,:,jans+ 9),3)*12/nyears)*md.materials.rho_freshwater/md.materials.rho_ice)'/1000;
smb.nov = ((sum(smbmonth(:,:,jans+10),3)*12/nyears)*md.materials.rho_freshwater/md.materials.rho_ice)'/1000;
smb.dec = ((sum(smbmonth(:,:,jans+11),3)*12/nyears)*md.materials.rho_freshwater/md.materials.rho_ice)'/1000;

% Interpolate onto model mesh
for fn=fieldnames(smb)'
smb_mesh.(fn{1})=InterpFromGridToMesh(x,y,smb.(fn{1}),md.mesh.x,md.mesh.y,0);
end

% Load and interpolate the annual smb for 1960-89
smb_6089 = ((sum(smbmonth,3)/nyears)*md.materials.rho_freshwater/md.materials.rho_ice)'/1000;
smb_6089 = InterpFromGridToMesh(x,y,smb_6089,md.mesh.x,md.mesh.y,0);

% Find difference between monthly maps and the mean yearly SMB
%    (smb_6089)
for fn=fieldnames(smb)'
  diff.(fn{1}) = smb_mesh.(fn{1}) - smb_6089;
end

clear smbmonth smb jans smb_mesh
%%}}}

% SMB: find seasonality (for 1960-1979 reference period) %%{{{
% NOTE: This is used to calculate the modeled dynamic thickness change anomaly in a way that is
%       consistent with SERAC thickness change observations, which also use 1960-1979 as a reference period
nyears = 20;
nmonths = nyears*12;
% monthly data from 20 years from Jan 1960 - Dec 1979:
smbmonth   = ncread(ncsmb,'SMB_rec',[1 1 25],[Inf Inf nmonths]); % monthly smb in mm WE / month

% Load and interpolate the annual smb for 1960-79
smb_6079 = ((sum(smbmonth,3)/nyears)*md.materials.rho_freshwater/md.materials.rho_ice)'/1000;
smb_6079 = InterpFromGridToMesh(x,y,smb_6079,md.mesh.x,md.mesh.y,0);

clear smbmonth smb ncsmb x y jans smb_mesh
%%}}}

% Loop through ensemble members to calculate SMB anomaly for each one %%{{{
for i = [1:9 26:153] % skip A0009 thru A0024

	id = lhs.IDs(i);

   fprintf('processing %s\n', id{1});

	% smb   
   % smb1: shift in mean smb
   %       up to ~+/- 127.5 Gt averaged over entire ice sheet (in m/yr ice)
   %       This is +30%/-30% of 1960-89 mean (or approximately 1 sigma)
   % 1. find ice sheet area:
   mdarea = sum(GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y));
   % 2. Read in LHS scalar and convert to +/-0.3
   % [0 1] --> [-30% +30%] (1 sigma of 1960-89 period)
   % [0 1] --> [-0.3 0.3]
   lhsv3 = -0.3 + lhs.smb1(i)*0.6;
   % 3. find ensemble value for smb1 (in Gt)
   mean6089 = findtotalSMB(md,smb_6089);
   smb1_Gt = lhsv3*mean6089;
   % 4. find m/yr anomaly (convert Gt to tons W.E. to m^3 of ice, then divide by area in m^2) 
   smb1anomaly = ((smb1_Gt*1e9)*md.materials.rho_freshwater/md.materials.rho_ice)/mdarea;
   % MIGHT WANT TO WEIGHT THIS DEPENDING ON REGION (+127.5 Gt distributed
   % over the whole ice sheet is +0.08 m/yr everywhere)
   % 5. add anomaly to the default annual smb map (2001-2015 average)
   pert_smb1 = smb_0115 + smb1anomaly;

   % smb2: varying strength of seasonal cycle
   %       from mean seasonal cycle of the 1960-1989 period
   %       mean cycle = 1
   %       varies between 0 (no cycle) and 2 (double amplitude)
   % 1. Find monthly maps of SMB (averaged for each month for 2001-15
   %    period --> done, see smb_mesh.XXX
   % 2. Find difference between monthly maps and the mean yearly SMB
   %    (smb_mesh.year) --> done, see diff
   % 3) Apply factor to these monthly difference maps (between 0 and 2)
   lhsv4 = lhs.smb2(i)*2;
   for fn=fieldnames(diff)'
    	pert.(fn{1}) = diff.(fn{1})*lhsv4;
   end
   % 4) Add back onto yearly mean map (what ever it is after the smb1
   %    anomaly has been added).
   for fn=fieldnames(diff)'
    	pert_smb2.(fn{1}) = pert.(fn{1}) + pert_smb1;
   end

   pert_smb = zeros(md.mesh.numberofvertices+1,12);
   for j = 1:12
    	fn = fieldnames(diff)';
    	mnth = fn(j);
    	pert_smb(1:(end-1),j) = [pert_smb2.(mnth{1})];
    	pert_smb(end,j) = ((j-1)/12);
   end

   pert_smb_trapz = trapz(pert_smb(end,:), pert_smb(1:end-1,:), 2); % m ice eq. over the year
   anom_smb = pert_smb_trapz - smb_6079;

   % output netcdf
   output_netcdf_filename = sprintf('/Volumes/LaCie/GrIS_Calibrated_SLR/committedSLR/netcdfs/gris.smb_anom6079.%s.nc', lhs.IDs{i});
	%output_netcdf_filename = sprintf('../netcdfs/gris.smb_anom.%s.test.nc', lhs.IDs{i});
	%if exist(output_netcdf_filename, 'file')
	%   fprintf([output_netcdf_filename ' already exists ... skipping!\n'])
	%   continue
	%end

	xgrid =  -720000:1000: 960000;
	ygrid = -3450000:1000:-570000;
	anom_smb_grid = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y, anom_smb, xgrid, ygrid, NaN);

	nccreate(output_netcdf_filename, 'smb_anom', 'Dimensions', {'y', length(ygrid), 'x', length(xgrid)});
	ncwrite( output_netcdf_filename, 'smb_anom', anom_smb_grid);
	nccreate(output_netcdf_filename, 'x', 'Dimensions', {'x', length(xgrid)});
	ncwrite( output_netcdf_filename, 'x', xgrid);
	nccreate(output_netcdf_filename, 'y', 'Dimensions', {'y', length(ygrid)});
	ncwrite( output_netcdf_filename, 'y', ygrid);

	fprintf('\n')

end %%}}}

