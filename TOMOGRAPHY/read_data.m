function [ data ] = read_data( DataFile,StaFile,par )

%=================== input rays =======================================
%
fid = fopen(DataFile, 'r');
A = textscan(fid, '%7.4f  %7.3f  %7.2f  %s  %4.0f  %8.4f  %9.4f  %7.3f  %6.2f  %6.4f  %6.4f'); 
fclose(fid);

ray.pd    = A{1};
ray.p     = r2d(ray.pd)/6371;
ray.gcarc = A{2};
ray.seaz  = A{3};
ray.sta   = A{4};
ray.orid  = A{5};
ray.elat  = A{6};
ray.elon  = A{7};
ray.edep  = A{8};
ray.d     = A{9};
ray.sd    = A{10};
ray.cf    = A{11};
clear A

ray.nrays = length(ray.p);

data.ray=ray;

%======================= save event info ==============================

[evt.orid, I, J] = unique(ray.orid);
evt.lat   = ray.elat(I);
evt.lon   = ray.elon(I);
evt.depth = ray.edep(I);
evt.nevts = length(evt.orid);

data.evt = evt;

%====================== input station =================================
sta_names = unique( ray.sta );

fid = fopen(StaFile,'r');
B = textscan(fid,'%s %f %f %f %s %s %f %f %f %s %f %f %f','delimiter',','); 
fclose(fid);

% stas{is}, slats(is), slons(is), selevs(is), statype{is}, IP,trm,...
%         xrdg,ocage,plate,yr,moh);

stns.sta = strtok(B{1},' ');  
stns.lat = B{2};  
stns.lon = B{3};  
stns.elv = B{4};
stns.typ = strtok(B{5},' ');  
stns.ip  = B{6}; 
stns.trm = B{7};  
stns.Xrd = B{8};
stns.age = B{9};
stns.plt = strtok(B{10},' ');  
stns.yr  = B{11};
stns.moh = B{12};
stns.sed = B{13};
clear B


stn.lat = zeros(length(sta_names),1);
stn.lon = zeros(length(sta_names),1);
stn.elv = zeros(length(sta_names),1);
stn.typ = cell(length(sta_names),1);
stn.ip  = cell(length(sta_names),1);
stn.trm = zeros(length(sta_names),1);
stn.Xrd = zeros(length(sta_names),1);
stn.age = zeros(length(sta_names),1);
stn.plt = cell(length(sta_names),1);
stn.yr  = zeros(length(sta_names),1);
stn.moh = zeros(length(sta_names),1);
stn.sed = zeros(length(sta_names),1);

% ------- associate data ------------------------------------------------
for ii = 1:length(sta_names)
    ix=find(strcmp(sta_names{ii}, stns.sta ));
    stn.lat(ii) = stns.lat(ix);
    stn.lon(ii) = stns.lon(ix);
    stn.elv(ii) = stns.elv(ix);
    stn.typ(ii) = stns.typ(ix);
    stn.ip(ii)  = stns.ip(ix);
    stn.trm(ii) = stns.trm(ix);
    stn.Xrd(ii) = stns.Xrd(ix);
    stn.age(ii) = stns.age(ix);
    stn.plt(ii) = stns.plt(ix);
    stn.yr(ii)  = stns.yr(ix);
    stn.moh(ii) = stns.moh(ix);
    stn.sed(ii) = stns.sed(ix);
    stn.num(ii) = ii;
end

stn.sta = sta_names;
stn.nstas = length(sta_names);

ray.sta_num = zeros(length(ray.d),1);
for jj = 1:length(stn.sta)
    ix=find( strcmp(stn.sta{jj},ray.sta) );
    ray.sta_num(ix) = jj;
    ray.stalat(ix,1) = stn.lat(jj);
    ray.stalon(ix,1) = stn.lon(jj);
end


% data structure
data.ray = ray;
data.stn = stn;

% add the estimate of the total ray lengths
data = trdist_est(data, par);



function [ data ] = trdist_est( data, par )
%crude estimate works well enough
% these rayp-distance values are from taup, using IASP91

if par.PS==1
    ipd   = [ 0     4.548   4.638   5.403  6.148  6.877  7.602  8.304  8.615  8.835  9.099  10.898 13.627 13.700];
elseif par.PS==2
    ipd   = [ 0     8.658   9.199   10.519 11.724 12.868 13.963 14.958 15.368 15.669 15.961 20.047 24.314 24.560];
end
idist     = [ 13000 10563.5 10007.5 8895.6 7783.6 6672   5560   4448   3892   3336   2780   2224   1668   1112];

pd = data.ray.pd;
data.ray.trdist = interp1(ipd,idist,pd);

end

end