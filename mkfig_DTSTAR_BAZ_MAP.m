% cycle through events and calculate differential travel times for all stations
clear all
close all
cd ~/Documents/MATLAB/CASC_atten/
addpath('matguts')

%% parameters
ifsave = 1;
method = 'comb';% 'comb' or 'specR'

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

cmap = parula;
% tstlim = 2.*[-1.5 0.5];

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
tstlim = ip*[-0.8 0.8];

% LOAD RESULTS
if strcmp(method,'specR')
    load([resdir,'all_dtstar_',phase,'_',component,'.mat']);
elseif strcmp(method,'comb')
    load([resdir,'all_dtstar',method,'_',phase,'_',component,'.mat']);
    all_dtstar = all_dtstar_comb;
else
    error('Need to specify method')
end

%% parse the stations and their data - collate repeats and remove allnan stas
[ stas,all_dtstar,slats,slons,selevs ] = results_PARSE_STATIONS( stas_all,all_dtstar,slats_all,slons_all,selevs_all );
nstas = length(stas);

%% ------------------ MAP WITH DT FOR THIS EVENT ------------------

figure(32), clf, hold on
mkfig_CascMAP
lonlims = [-132.1 -120];
set(gcf,'position',[200 300 plotsize/plot_size_ratio(lonlims,latlims) plotsize])
set(gca,'xlim',lonlims)

% nnanstas = ~isnan(nanmean(all_dT,2));
Nobs_sta = sum(~isnan(all_dtstar),2);

% loop over stations
for is = 1:nstas
    nnan = find(~isnan(all_dtstar(is,:)));
    for ie = 1:length(nnan)
        baz = azimuth(slats(is),slons(is),elats(nnan(ie)),elons(nnan(ie)));
        p1la = slats(is);
        p1lo = slons(is);
        [p2la,p2lo] = reckon(p1la,p1lo,scale,baz);
        % plot line toward baz, coloured by tstar
        plot([p1lo;p2lo],[p1la;p2la],'LineWidth',2,...
            'color',colour_get(all_dtstar(is,nnan(ie)),tstlim(2),tstlim(1),cmap))
    end % loop on nnan evts
end % loop on stas

%% colour bar
colormap(cmap), caxis(tstlim)
tkvl = unique(round_level([tstlim(1):0.5:tstlim(2)],0.5));
cbar_custom(gca, 'location',[-131.4 -131. 39.2 42.5],'tickside','right',...
    'lims',tstlim,'tickvals',tkvl,'cmap',cmap,...
    'FontSize',12,'FontWeight','bold',...
	'title',sprintf('$\\Delta t^*_%s$ \\,(s)',phase),'interpreter','latex');

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
title(sprintf('$\\Delta t^*$ for $%s$-waves (%s component) %s',phase,component,method),...
      'FontSize',18,'FontWeight','bold','Interpreter','latex')
set(gca,'FontSize',14,'LineWidth',2.5,'box','on')

% save
if ifsave
save2pdf(32,sprintf('dtstar_map_baz_%s_%s_%s',method,phase,component),'figs');
end

end