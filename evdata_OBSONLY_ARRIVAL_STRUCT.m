% We want to do the differential TT and differential TSTAR measurements on
% OBS stations only, so here is a script to replicate eqar arrival
% structures and keep only the OBS stations. These structures can then be
% straightforwardly fed into calc_COMBSPECTRA or calc_DIFFERENTIAL_TT

clear all
addpath('matguts')

%% Parameters
phase = 'P';
component = 'Z';

%% Directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% DATA DIRECTORY (top level)
datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash

%% GET STATION + EVENT DATA
[ nstas,stas,slats,slons,selevs,ondate,offdate,staname,statype ] = db_stadata( dbdir,dbnam );
[ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam ) ;
  
 %% Loop over events
 
for ie = 1:norids % % loop on orids
    fprintf('\n Orid %.0f %s \n\n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
    % name files and directories
    evdir       = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];

    INdatinfofile = [datadir,evdir,'_datinfo_',phase];
    INarfile      = [datadir,evdir,'_EQAR_',phase,'_',component];

    OUTdatinfofile = [datadir,evdir,'_datinfo_OBS_',phase];
    OUTarfile      = [datadir,evdir,'_EQAR_OBS_',phase,'_',component];
    
    % check files exist
    if ~exist([INdatinfofile,'.mat'],'file'), fprintf('No station mat files for this event\n');continue, end
    if ~exist([INarfile,'.mat'],'file'), fprintf('No arfile for this event + phase\n');continue, end
    
    % load files
    load(INdatinfofile) % loads datinfo stucture
    load(INarfile)      % loads eqar structure

    %% Get station names - parse OBS
    if ~isequal(length({datinfo.sta}'),length({eqar.sta}'))
        error('datinfo and eqar files do not have same stations')
    end
    
    ar_stas = {datinfo.sta}';
    [~,~,arinds] = intersect(ar_stas,stas,'stable');
    gdstas = strcmp('OBS',statype(arinds));

    if ~any(gdstas)
        fprintf('No OBS stations for this EQ\n')
        continue
    else
        datinfo = datinfo(gdstas);
        eqar = eqar(gdstas);
        save(OUTdatinfofile,'datinfo')
        save(OUTarfile,'eqar')
        fprintf(' only OBS arrival files saved\n')        
    end

end

 
 