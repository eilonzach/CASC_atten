% cycle through events and calculate differential travel times for all stations
clear all
close all
cd ~/Documents/MATLAB/CASC_atten/
addpath('matguts')

%% parameters
ifsave = 1;
method = 'xcorr';% 'comb' or 'xcorr' 

plotsize = 800;
scale = 0.35; % length of lines
tdftick = 1;
phases = {'P','S'};
components = {'Z','T'}; %'Z', 'R', or 'T'

keyloc = [-131.1,49.9];
keysiz = [1.4,1.9];

%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% DATA DIRECTORY (top level)
datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash
% RESULTS DIRECTORY
resdir = '~/Documents/MATLAB/CASC_atten/results/'; % needs final slash
addpath('~/Documents/MATLAB/seizmo-master/cmap/');

cmap = blue2red; 
close all
%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

% GET EVENTS+STATIONS DATA
[ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam );
[ nstas_all,stas_all,slats_all,slons_all,selevs_all ] = db_stadata( dbdir,dbnam );

% PARSE RESULTS
% results_parse

%% Get results
for ip = 1:length(phases)
phase = phases{ip};
component = components{ip};
tdflim = ip*[-1.1 1.1];

% LOAD RESULTS
if strcmp(method,'xcorr')
    load([resdir,'all_dT_',phase,'_',component,'.mat']);
elseif strcmp(method,'comb')
    load([resdir,'all_dT',method,'_',phase,'_',component,'.mat']);
    all_dT = all_dT_comb;
else
    error('Need to specify method')
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
all_dT = all_dT + topo_cor_S*ones(1,size(all_dT,2));

%% ------------------ MAP WITH DT FOR THIS EVENT ------------------

figure(32), clf, hold on
mkfig_CascMAP
lonlims = [-132.1 -120];
set(gcf,'position',[200 300 plotsize/plot_size_ratio(lonlims,latlims) plotsize])
set(gca,'xlim',lonlims)

% nnanstas = ~isnan(nanmean(all_dT,2));
Nobs_sta = sum(~isnan(all_dT),2);

% loop over stations
for is = 1:nstas
    nnan = find(~isnan(all_dT(is,:)));
    for ie = 1:length(nnan)
        baz = azimuth(slats(is),slons(is),elats(nnan(ie)),elons(nnan(ie)));
        p1la = slats(is);
        p1lo = slons(is);
        [p2la,p2lo] = reckon(p1la,p1lo,scale,baz);
        % plot line toward baz, coloured by arrival time
        plot([p1lo;p2lo],[p1la;p2la],'LineWidth',2,...
            'color',colour_get(all_dT(is,nnan(ie)),tdflim(2),tdflim(1),cmap))
    end % loop on nnan evts
end % loop on stas

%% colour bar
colormap(cmap), caxis(tdflim)
tkvl = unique(round_level([tdflim(1):ip*0.5:tdflim(2)],0.5));
cbar_custom(gca, 'location',[-131.4 -131.0 39.2 42.5],'tickside','right',...
    'lims',tdflim,'tickvals',tkvl,'cmap',cmap,...
    'FontSize',12,'FontWeight','bold',...
	'title',sprintf('$\\delta T_%s$ \\,(s)',phase),'interpreter','latex');

% %% key
% kle = keyloc(1) - 0.5*keysiz(1); kri = keyloc(1) + 0.5*keysiz(1);
% kbo = keyloc(2) - 0.5*keysiz(2); kto = keyloc(2) + 0.5*keysiz(2);
% kw = keysiz(1); kh = keysiz(2);
% 
% hold on
% patch([kle kle kri kri kle],[kbo kto kto kbo kbo],...
%       'w','LineWidth',2)
% text(keyloc(1),kto-0.1*kh,'Nevts',...
%     'fontsize',15,'fontweight','bold','verticalalignment','top','horizontalalignment','center');
% 
% scatter(kle + 0.3*kw,kbo + 0.61*kh,scale*sqrt(1),  'k','filled','MarkerEdgeColor','k');
% scatter(kle + 0.3*kw,kbo + 0.435*kh,scale*sqrt(10), 'k','filled','MarkerEdgeColor','k');
% scatter(kle + 0.3*kw,kbo + 0.18*kh,scale*sqrt(50),'k','filled','MarkerEdgeColor','k');
% 
% text(kle + 0.62*kw,kbo + 0.61*kh,'1', 'fontsize',12,'fontweight','bold','verticalalignment','middle');
% text(kle + 0.62*kw,kbo + 0.435*kh,'10','fontsize',12,'fontweight','bold','verticalalignment','middle');
% text(kle + 0.62*kw,kbo + 0.18*kh,'50','fontsize',12,'fontweight','bold','verticalalignment','middle');

%% title
title(sprintf('$\\delta T$ for $%s$-waves (%s component)',phase,component),...
      'FontSize',18,'FontWeight','bold','Interpreter','latex')
set(gca,'FontSize',14,'LineWidth',2.5,'box','on')

% save
if ifsave
save2pdf(32,sprintf('dT_map_baz_%s_%s_%s',method,phase,component),'figs');
end

end % loop on phases
