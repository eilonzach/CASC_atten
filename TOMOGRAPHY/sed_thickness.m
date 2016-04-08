function [ sta_sed ] = sed_thickness( slats,slons )
%SED_THICKNESS Summary of this function goes here
%   Function to output the sediment thickness for OBS stations, based on
%   interpolated points on a sediment grid file

obs_sed_file = '/Users/zeilon/Work/CASCADIA/plots/seds/obs_seds.txt';
fid = fopen(obs_sed_file,'r');
A = textscan(fid,'%f %f %f','Headerlines',1);
fclose(fid);
sed_lons = A{1};
sed_lats = A{2};
sed_hs = A{3};

F = scatteredInterpolant(sed_lons,sed_lats,sed_hs);

sta_sed = F(slons,slats);


end

