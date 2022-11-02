function [ totalsmb ] = findtotalSMB( md, smb )
%FINDTOTALSMB Summary of this function goes here
%   Detailed explanation goes here
    
    
    index = md.mesh.elements;
    x = md.mesh.x;
    y = md.mesh.y;
    
    % 1. get some parameters
    rho_ice   = md.materials.rho_ice;
    rho_freshwater = md.materials.rho_freshwater;

    % 2. compute averages across each element
    smbmean  = mean(smb(index),2);
    
    % 3. get areas of all triangles (elements)
    areas = GetAreas(index,x,y);
    
    % 4. Compute volume of smb in each element
    V = areas.*smbmean;
    
    % 5. take out the ones that are outside of levelset
    %pos = find(md.mask.ice_levelset(index),[],2)>0));
    %V(pos) = 0;

    % sum individual contributions
    V = sum(V);

    % convert from m^3/yr ice to Gt/yr WE
    totalsmb = (V*rho_ice/rho_freshwater)/1e9;
    
    
end

