%% parse results into large nstas*nevts structures
clear all
% if running alone, establish these, otherwise use previous values
phase = 'P';
component = 'Z';

%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% DATA DIRECTORY (top level)
datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash
% RESULTS DIRECTORY
resdir = '~/Documents/MATLAB/CASC_atten/TOMOGRAPHY/data/'; % needs final slash
addpath('/Users/zeilon/Documents/MATLAB/CASC_atten/matguts')

%% out files
dT_file = sprintf('data_dT_%s%s.dat',phase,component);
dtstar_file = sprintf('data_dtstar_%s%s.dat',phase,component);
sta_file = 'stations.dat';

%% conditions
mag_min = 6.25; % skip events if below this mag
acor_min = 0.8; % skip events if xcor max below this acor
snr_min = 10; % skip result if data below this snr

%% parms
dTsd = 0.05; % standard minimum standard deviation. Divide this by acor
dtssd = 0.1; % standard minimum standard deviation. Divide this by acor


%% db data
[ nstas,stas,slats,slons,selevs,ondate,offdate,staname,statype ] = db_stadata( dbdir,dbnam );
[ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam );
sages = jdf_crust_age(slats,slons,'linear');

%% =================================================================== %%
%% ============================  STAS  =============================== %%
%% =================================================================== %%

% fields: sta  slat  slon  selev  obs/land  OBS_IP  ifTRM  xrdg  oceanic_age  jdforgor  yr  moho
fid = fopen([resdir sta_file],'w');
for is = 1:nstas
    % do things    
    
    [ dist_deg,xrdg,plate ] = dist_to_rdg( slats(is),slons(is) );
    
    if strcmp(statype(is),'OBS')
    
        if plate==1; jdf_or_gor='jdf';elseif plate==0, jdf_or_gor='gor'; end
        
        [IP,note] = which_OBS(stas{is}); 
        if strcmp('TRM',note), trm = 1; else trm=0; end
        yr=0;
        switch stas{is}(end) 
            case 'A', yr=1;
            case 'B', yr=2;
            case 'C', yr=3;
            case 'D', yr=4;
        end
        
        moh = 6.0;
        
        sed = sed_thickness(slats(is),slons(is));
        if isnan(sed), sed=-999; end

        ocage = sages(is);
        
    else % land defaults
        
        IP = []; trm = 0; yr = 0; moh = 0; sed = 0; jdf_or_gor=[]; ocage = 0;
    
    end
    
    % sort out annoying double-named stations
    sta = which_CASC_STA( stas{is},slats(is),slons(is) );
    
    fprintf(fid,'%-6s, %8.4f, %9.4f, %7.1f, %-4s, %4s, %1.0f, %6.1f, %5.2f, %3s, %1.0f, %5.2f, %5.2f\n',...
        sta, slats(is), slons(is), selevs(is), statype{is}, IP,trm,...
        xrdg,ocage,jdf_or_gor,yr,moh,sed);
end
fclose(fid) ;   

fprintf('Not doing data\n')

%% =================================================================== %%
%% ===========================  DATA  ============================== %%
%% =================================================================== %%
dT_data = struct('rayp',[],'gcarc',[],'seaz',[],'sta',{},...
                  'orid',[],'elat',[],'elon',[],'edep',[],...
                  'dat',[],'sd',[],'cfreq',[]);
dtstar_data = struct('rayp',[],'gcarc',[],'seaz',[],'sta',{},...
                  'orid',[],'elat',[],'elon',[],'edep',[],...
                  'dat',[],'sd',[],'cfreq',[]);
              

count_dT = 0;
count_dts = 0;
              
for ie = 1:402 % 44:norids % loop on orids
    fprintf('Orid %.0f\n',ie)
    if  mags(ie)<mag_min, continue, end
    
    %% grab data!
    % name files and directories
    evdir       = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
    datinfofile = [datadir,evdir,'_datinfo_',phase];
    arfile      = [datadir,evdir,'_EQAR_',phase,'_',component];
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
    if length([eqar.dT])~=length(eqar),for is = 1:length(eqar),if isempty(eqar(is).dT),eqar(is).dT=nan;end,end,end
    dTgind = find(~isnan([eqar.dT]));
    dTgind = dTgind([eqar(dTgind).acor] >= acor_min);
    mean_dT = mean([eqar(dTgind).dT]);
    for is = 1:length(dTgind)
        ind = dTgind(is);
        count_dT = count_dT+1;
        
        dT_data(count_dT,1) ...
         = struct('rayp',eqar(ind).rayp,...
                  'gcarc',eqar(ind).gcarc,...
                  'seaz',eqar(ind).seaz,...
                  'sta',which_CASC_STA(eqar(ind).sta,eqar(ind).slat,eqar(ind).slon),...
                  'orid',orids(ie),'elat',elats(ie),'elon',elons(ie),'edep',edeps(ie),...
                  'dat',eqar(ind).dT - mean_dT,...
                  'sd',dTsd./eqar(ind).acor,...
                  'cfreq',mean([eqar(ind).par_dT.fhi,eqar(ind).par_dT.flo]));
    end
    
    %% parse dtstar_COMB data
    if ~isfield(eqar,'dtstar_comb'), fprintf('   NEED TO CALC DTSTAR w COMB\n',phase), continue, end
    if length([eqar.dtstar_comb])~=length(eqar), error('Not enough dtstars somehow'), end
    dtsgind = find(~isnan([eqar.dtstar_comb]));
    mean_dtstar = mean([eqar(dtsgind).dtstar_comb]);
    for is = 1:length(dtsgind)
        ind = dtsgind(is);
        count_dts = count_dts+1;
        
        dtstar_data(count_dts,1) ...
         = struct('rayp',eqar(ind).rayp,...
                  'gcarc',eqar(ind).gcarc,...
                  'seaz',eqar(ind).seaz,...
                  'sta',which_CASC_STA(eqar(ind).sta,eqar(ind).slat,eqar(ind).slon),...
                  'orid',orids(ie),'elat',elats(ie),'elon',elons(ie),'edep',edeps(ie),...
                  'dat',eqar(ind).dtstar_comb - mean_dtstar,...
                  'sd',dtssd,...
                  'cfreq',mean([1./eqar(ind).par_dtstar_comb.comb.Tmin,1./eqar(ind).par_dtstar_comb.comb.Tmax]));
    end

end % loop on orids

fprintf('%.0f dtstar measurements for %s on %s comp\n',count_dts,phase,component)
fprintf('%.0f dT measurements for %s on %s comp\n',count_dT,phase,component)

%% =================================================================== %%
%% =============================  dT  ================================ %%
%% =================================================================== %%

% fields: rayp  gcarc  seaz  sta  orid  elat  elon  edep  dT  sd  c_freq
fid = fopen([resdir dT_file],'w');
for ii = 1:count_dT
fprintf(fid,'%7.4f  %7.3f  %7.2f  %6s  %4.0f  %8.4f  %9.4f  %7.3f  %6.2f  %6.4f  %6.4f\n',...
    dT_data(ii).rayp,dT_data(ii).gcarc,dT_data(ii).seaz,dT_data(ii).sta,...
    dT_data(ii).orid,dT_data(ii).elat,dT_data(ii).elon,dT_data(ii).edep,...
    dT_data(ii).dat,dT_data(ii).sd,dT_data(ii).cfreq);
end
fclose(fid);
   


%% =================================================================== %%
%% ===========================  dtstar  ============================== %%
%% =================================================================== %%

% fields: rayp  gcarc  seaz  sta  orid  elat  elon  edep  dtstar  sd  c_freq
fid = fopen([resdir dtstar_file],'w');
for ii = 1:count_dts
fprintf(fid,'%7.4f  %7.3f  %7.2f  %6s  %4.0f  %8.4f  %9.4f  %7.3f  %6.2f  %6.4f  %6.4f\n',...
    dtstar_data(ii).rayp,dtstar_data(ii).gcarc,dtstar_data(ii).seaz,dtstar_data(ii).sta,...
    dtstar_data(ii).orid,dtstar_data(ii).elat,dtstar_data(ii).elon,dtstar_data(ii).edep,...
    dtstar_data(ii).dat,dtstar_data(ii).sd,dtstar_data(ii).cfreq);
end
fclose(fid);



       

    
    
    
    