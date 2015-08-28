% cycle through events and calculate spectral ratios for all stations
clear all
addpath('matguts')

%% parameters

%% directories 
% RESULTS DIRECTORY
rdir = '~/Documents/MATLAB/CASC_atten/results_SPEC_RATIOS/'; % needs final slash
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascattendb';
% DATA DIRECTORY (top level)
datadir = '/VolumeS/DATA/CASCADIA/DATA/'; % needs final slash

%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

% GET EVENTS DATA
db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,elats,elons,edeps,evtimes,mags] = dbgetv(dbor,'orid','lat','lon','depth','time','ms');
norids = dbnrecs(dbor);
dbclose(db);


for ie = 1:100 % 1:norids % loop on orids
    fprintf('\n Orid %.0f %s \n\n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
    evdir = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
    if ~exist([datadir,evdir,'_datinfo.mat'],'file'), fprintf('No station mat files for this event\n');continue, end
    load([datadir,evdir,'_datinfo.mat'])
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end
    
    return
    % RESULTS STRUCTURE
%     eq = struct('
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end