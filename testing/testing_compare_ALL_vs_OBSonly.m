% cycle through events and calculate differential travel times for all stations
clear all
close all
cd ~/Documents/MATLAB/CASC_atten/
addpath('matguts')

%% parameters
ifsave = 1;

phases = {'P','S'};
components = {'Z','T'}; %'Z', 'R', or 'T'
method_dtstar = 'comb';% 'comb' or 'specR'

plotsize = 800;
scale = 100; % length of lines

%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% DATA DIRECTORY (top level)
datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash
% RESULTS DIRECTORY
resdir = '~/Documents/MATLAB/CASC_atten/results/'; % needs final slash
addpath('~/Documents/MATLAB/seizmo-master/cmap/');


% GET EVENTS+STATIONS DATA
[ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam );
[ nstas_all,stas_all,slats_all,slons_all,selevs_all,~,~,~,statype ] = db_stadata( dbdir,dbnam );

% PARSE RESULTS
% results_parse

%% Get results
for ip = 1:length(phases)
phase = phases{ip};
component = components{ip};
tstlim = ip*[-0.8 0.8];

%% DTSTAR
% LOAD DTSTAR RESULTS
if strcmp(method_dtstar,'specR')
    load([resdir,'all_dtstar_OBS_',phase,'_',component,'.mat']);
    obs_dtstar = all_dtstar;
    load([resdir,'all_dtstar_',phase,'_',component,'.mat']);
    full_dtstar = all_dtstar;
elseif strcmp(method_dtstar,'comb')
    load([resdir,'all_dtstar',method_dtstar,'_OBS_',phase,'_',component,'.mat']);
    obs_dtstar = all_dtstar_comb;
    load([resdir,'all_dtstar',method_dtstar,'_',phase,'_',component,'.mat']);
    full_dtstar = all_dtstar_comb;
else
    error('Need to specify method')
end

% limit to OBS
nobs = sum(strcmp(statype,'OBS'));
obs_dtstar = obs_dtstar(strcmp(statype,'OBS'),:);
full_dtstar = full_dtstar(strcmp(statype,'OBS'),:);
% take off event mean
obs_dtstar = obs_dtstar - ones(nobs,1)*nanmean(obs_dtstar,1);
full_dtstar = full_dtstar - ones(nobs,1)*nanmean(full_dtstar,1);



%% DT
% LOAD DT RESULTS
load([resdir,'all_dT_OBS_',phase,'_',component,'.mat']);
obs_dT = all_dT;
load([resdir,'all_dT_',phase,'_',component,'.mat']);
full_dT = all_dT;

% limit to OBS
obs_dT = obs_dT(strcmp(statype,'OBS'),:);
full_dT = full_dT(strcmp(statype,'OBS'),:);
% take off event mean
obs_dT = obs_dT - ones(nobs,1)*nanmean(obs_dT,1);
full_dT = full_dT - ones(nobs,1)*nanmean(full_dT,1);


%% plots
figure(23), hold on
scatter(obs_dtstar(:),full_dtstar(:),30,'filled')
plot([-1 1],[-1 1],'--k')
title('DTSTAR','FontSize',19)
xlabel('OBS only','FontSize',16)
ylabel('ALLSTAS','FontSize',16)



figure(24), hold on
scatter(obs_dT(:),full_dT(:),30,'filled')
plot([-1 1],[-1 1],'--k')
title('DT','FontSize',19)
xlabel('OBS only','FontSize',16)
ylabel('ALLSTAS','FontSize',16)



return
end
