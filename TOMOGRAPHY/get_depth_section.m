function [dV,anis,hq,slt,sln,ss,Dxy,smbv,smba,Xxy] = get_depth_section(pt1,pt2,plot_model)
% inputs:
%   pt1     = [lon1,lat1];
%   pt2     = [lon2,lat2];
% outputs:
%   val      = matrix of vals within slice
%   hq      = matrix of hit quality within slice
%   slt     = vector of latitude along section
%   sln     = vector of longitude along section
%   ss      = vector of distance (in degrees) along section
%   Dxy     = total distance (in degrees) along section
%   smbv    = matrix of velocity semblance within slice
%   smba    = matrix of anisotropy semblance within slice
%   Xxy    = total distance (in km) along section


pmlt = plot_model.lt(:,:,1);
pmln = plot_model.ln(:,:,1);

%% Size and location of section
lat1 = pt1(2); lon1 = pt1(1);
lat2 = pt2(2); lon2 = pt2(1);
Dxy = distance(lat1,lon1,lat2,lon2); %in degrees

%% Get coords of section
nxy = 20*round(Dxy);
sln = linspace(lon1,lon2,nxy)'; % lon points along profile
slt = linspace(lat1,lat2,nxy)'; % lat points along profile
ss  = linspace(0,Dxy,nxy)'; % distance along profile

nz = length(unique(plot_model.z));

%% get values of model along section
val = zeros(nz,nxy);
hq = zeros(nz,nxy);
for iz = 1:nz
    val(iz,:)   = griddata(pmln,pmlt,plot_model.val(:,:,iz),sln,slt);
    hq(iz,:)   = griddata(pmln,pmlt,plot_model.hq(:,:,iz),sln,slt);
end

%% get the semblance values too
smbv = zeros(nz,nxy);
smba = zeros(nz,nxy);
for iz = 1:nz
    smbv(iz,:) = griddata(pmln,pmlt,plot_model.semb_v(:,:,iz),sln,slt);
    smba(iz,:) = griddata(pmln,pmlt,plot_model.semb_a(:,:,iz),sln,slt);
end

Xxy = distance_km(pt1(2),pt1(1),pt2(2),pt2(1));


