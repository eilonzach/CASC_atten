%% parse results into large nstas*nevts structures
clear all
% if running alone, establish these, otherwise use previous values
if ~exist('phase') || ~exist('component') 
phase = 'S';
component = 'T';
end

ifOBSonly = false;

%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% DATA DIRECTORY (top level)
datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash
% RESULTS DIRECTORY
resdir = '~/Documents/MATLAB/CASC_atten/results/'; % needs final slash

%% conditions
mag_min = 6.25; % skip events if below this mag
acor_min = 0.8; % skip events if xcor max below this acor
snr_min = 10; % skip result if data below this snr

%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

% GET EVENTS+STATIONS DATA
db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,elats,elons,edeps,evtimes,mags] = dbgetv(dbor,'orid','lat','lon','depth','time','ms');
norids = dbnrecs(dbor);
dbsi = dblookup_table(db,'site');
[stas,slats,slons,selevs] = dbgetv(dbsi,'sta','lat','lon','elev');
nstas = dbnrecs(dbsi);
dbclose(db)

obsstr = ''; if ifOBSonly, obsstr = 'OBS_'; end

all_dT     = nan(nstas,norids);
all_dtstar = nan(nstas,norids);
all_dT_comb = nan(nstas,norids);
all_dtstar_comb = nan(nstas,norids);

for ie = 1:402 % 44:norids % loop on orids
    fprintf('Orid %.0f\n',ie)
    if  mags(ie)<mag_min, continue, end
    %% grab data!
    % name files and directories
    evdir       = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
    datinfofile = [datadir,evdir,'_datinfo_',obsstr,phase];
    arfile      = [datadir,evdir,'_EQAR_',obsstr,phase,'_',component];
    % check files exist
    if ~exist([datinfofile,'.mat'],'file'), fprintf('No station mat files for this event\n');continue, end
    if ~exist([arfile,'.mat'],'file'), fprintf('No arfile for this event + phase\n');continue, end    
    % load info file
    load(datinfofile) % loads datinfo stucture
    % options to skip
    if isempty(datinfo), fprintf('No station mat files for this event\n'); continue, end

    % load data file
    load(arfile)      % loads eqar structure
    if all([eqar.pred_arrT]==0), fprintf('No %s arrivals for this event\n',phase), continue, end
    
%% parse dT data
    if ~isfield(eqar,'dT'), fprintf('   NEED TO XCOR\n',phase), continue, end
    for is = 1:length(eqar)
        if eqar(is).acor < acor_min, continue; end
        if strcmp(eqar(is).sta,'M08C'), eqar(is).slon = -124.895400; end
        if isempty(eqar(is).dT), continue, end
        [~,ind,~] = intersect(stas,datinfo(is).sta,'stable');
        all_dT(ind,ie) = eqar(is).dT;
    end
    
    %% parse dtstar data
    if ~isfield(eqar,'dtstar'), fprintf('   NEED TO CALC DTSTAR\n',phase), continue, end
    for is = 1:length(eqar)
        if eqar(is).snr_wf < snr_min, continue; end
        if strcmp(eqar(is).sta,'M08C'), eqar(is).slon = -124.895400; end
        if isempty(eqar(is).dtstar), continue, end
        [~,ind,~] = intersect(stas,datinfo(is).sta,'stable');
        all_dtstar(ind,ie) = eqar(is).dtstar;
    end
    
    %% parse dtstar_COMB data
    if ~isfield(eqar,'dtstar_comb'), fprintf('   NEED TO CALC DTSTAR w COMB\n',phase), continue, end
    for is = 1:length(eqar)
        if eqar(is).snr_wf < snr_min, continue; end
        if strcmp(eqar(is).sta,'M08C'), eqar(is).slon = -124.895400; end
        if isempty(eqar(is).dtstar_comb), continue, end
        if isnan(eqar(is).dtstar_comb), continue, end
        [~,ind,~] = intersect(stas,datinfo(is).sta,'stable');
        all_dT_comb(ind,ie) = eqar(is).dT_comb;
        all_dtstar_comb(ind,ie) = eqar(is).dtstar_comb;
    end
    
end % loop on orids

save([resdir,'all_dT_',obsstr,phase,'_',component],'all_dT')
save([resdir,'all_dtstar_',obsstr,phase,'_',component],'all_dtstar')
save([resdir,'all_dTcomb_',obsstr,phase,'_',component],'all_dT_comb')
save([resdir,'all_dtstarcomb_',obsstr,phase,'_',component],'all_dtstar_comb')


    
    
    
    