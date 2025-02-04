% cycle through events and calculate differential travel times for all stations
clear all
close all
cd ~/Documents/MATLAB/CASC_atten/
addpath('matguts')
addpath('~/Work/Codes/seizmo/cmap/')

%% parameters
ifsave = false;

method = 'xcorr';% 'comb' or 'xcorr' 

ifOBSinv = false;
ifOBSonly = true;

plotsize = 800;
scale = 100; % length of lines
tdftick = 1;
phases = {'P','S'};
components = {'Z','T'}; %'Z', 'R', or 'T'
ALP = []; % [] for the default (also 0), or the value of alpha

keyloc = [-131.1,49.9];
keysiz = [1.4,1.9];

%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% DATA DIRECTORY (top level)
datadir = '/Volumes/DATA_mini/CASCADIA/DATA/'; % needs final slash
% RESULTS DIRECTORY
resdir = '~/Documents/MATLAB/CASC_atten/results/'; % needs final slash
% FIGURES DIRECTORY
figdir = 'figs';

cmap = blue2red; 

close all
%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

obsloadstr = ''; 
if ifOBSinv, obsloadstr = 'OBS_'; obsstr = 'OBSinv_'; end

if isempty(ALP), alpstr = ''; else alpstr = sprintf('_%03.0f',ALP*100); end

% GET EVENTS+STATIONS DATA
[ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam );
[ nstas_all,stas_all,slats_all,slons_all,selevs_all,~,~,~,statypes_all ] = db_stadata( dbdir,dbnam );

% PARSE RESULTS
% results_parse

%% Get results
for ip = 1:length(phases)
phase = phases{ip};
component = components{ip};
tdflim = ip*[-1.1 1.1];

% LOAD RESULTS
if strcmp(method,'xcorr')
    load([resdir,'all_dT_',obsloadstr,phase,'_',component,alpstr,'.mat']); % alpstr has to be '' for this - just a check.
elseif strcmp(method,'comb')
    load([resdir,'all_dT',method,'_',obsloadstr,phase,'_',component,alpstr,'.mat']);
    all_dT = all_dT_comb;
else
    error('Need to specify method')
end

if ifOBSonly
    all_dT(~strcmp(statypes_all,'OBS'),:) = nan;
    obsstr = 'OBSonly_';
else 
    obsstr = '';
end
    

%% parse the stations and their data - collate repeats and remove allnan stas
[ stas,all_dT,slats,slons,selevs ] = results_PARSE_STATIONS( stas_all,all_dT,slats_all,slons_all,selevs_all );
nstas = length(stas);

%% topography correction
topo_cor_P = -selevs/6500; % add this to dT, it is positive if the elevation is negative (i.e. the station is missing travel time)
topo_cor_S = -selevs/3400; % add this to dT, it is positive if the elevation is negative (i.e. the station is missing travel time)
if strcmp(phase,'P'), topo_cor = topo_cor_P; elseif strcmp(phase,'S'), topo_cor = topo_cor_S; end
% sta_terms = sta_terms + topo_cor;
% sta_terms = sta_terms - mean(sta_terms);
all_dT = all_dT + topo_cor*ones(1,size(all_dT,2));

%% ------------------ MAP WITH DT FOR THIS EVENT ------------------

figure(31), clf, hold on
mkfig_CascMAP
lonlims = [-132.1 -120];
set(gcf,'position',[200 300 plotsize/plot_size_ratio(lonlims,latlims) plotsize])
set(gca,'xlim',lonlims)

% nnanstas = ~isnan(nanmean(all_dT,2));
Nobs = sum(~isnan(all_dT),2);
% compute station averages
[ sta_terms,evt_terms ] = lsq_sta_evt( all_dT,0.01 );



%% plot the station-averaged differential travel time
hp = scatter(slons,slats,scale*sqrt(Nobs),sta_terms,'filled','MarkerEdgeColor','k');

colormap(cmap), caxis(tdflim)

%% colour bar
tkvl = unique(round_level([tdflim(1):ip*0.5:tdflim(2)],0.5));
cbar_custom(gca, 'location',[-131.4 -131.0 39.2 42.5],'tickside','right',...
    'lims',tdflim,'tickvals',tkvl,'cmap',cmap,...
    'FontSize',12,'FontWeight','bold',...
	'title',sprintf('$\\delta T_%s$ \\,(s)',phase),'interpreter','latex');

%% key
kle = keyloc(1) - 0.5*keysiz(1); kri = keyloc(1) + 0.5*keysiz(1);
kbo = keyloc(2) - 0.5*keysiz(2); kto = keyloc(2) + 0.5*keysiz(2);
kw = keysiz(1); kh = keysiz(2);

hold on
patch([kle kle kri kri kle],[kbo kto kto kbo kbo],...
      'w','LineWidth',2)
text(keyloc(1),kto-0.1*kh,'Nevts',...
    'fontsize',15,'fontweight','bold','verticalalignment','top','horizontalalignment','center');

scatter(kle + 0.3*kw,kbo + 0.61*kh,scale*sqrt(1),  'k','filled','MarkerEdgeColor','k');
scatter(kle + 0.3*kw,kbo + 0.435*kh,scale*sqrt(10), 'k','filled','MarkerEdgeColor','k');
scatter(kle + 0.3*kw,kbo + 0.18*kh,scale*sqrt(50),'k','filled','MarkerEdgeColor','k');

text(kle + 0.62*kw,kbo + 0.61*kh,'1', 'fontsize',12,'fontweight','bold','verticalalignment','middle');
text(kle + 0.62*kw,kbo + 0.435*kh,'10','fontsize',12,'fontweight','bold','verticalalignment','middle');
text(kle + 0.62*kw,kbo + 0.18*kh,'50','fontsize',12,'fontweight','bold','verticalalignment','middle');

%% title
title(sprintf('%s $\\delta T$ for $%s$-waves (%s component)',strtok(obsstr,'_'),phase,component),...
      'FontSize',18,'FontWeight','bold','Interpreter','latex')
set(gca,'FontSize',14,'LineWidth',2.5,'box','on')

% save
if ifsave
save2pdf(31,sprintf('MAP_dT_staav_%s_%s%s_%s',method,obsstr,phase,component),figdir);

results = struct('stas',{stas},'dT',sta_terms,'slats',slats,'slons',slons,'selevs',selevs,'Nobs',Nobs);
resfile = sprintf('stav_dT%s_%s%s_%s',method,obsstr,phase,component);
save([resdir,resfile],'results')
end
pause
end % loop on phases
