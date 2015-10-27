% cycle through events and calculate differential travel times for all stations
clear all
close all
cd ~/Documents/MATLAB/CASC_atten/
addpath('matguts')

%% parameters
ifsave = 1;

scale = 100; % length of lines
tdftick = 1;
phases = {'P'};
component = 'Z'; %'Z', 'R', or 'T'

keyloc = [-131.1,49.9];
keysiz = [1.4,1.9];

%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascattendb';
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
[ nstas,stas,slats,slons,selevs ] = db_stadata( dbdir,dbnam );

% PARSE RESULTS
% results_parse

for ip = 1:length(phases)
phase = phases{ip};
tdflim = ip*1.5*[-1 1];

% LOAD RESULTS
load([resdir,'all_dT_',phase,'_',component,'.mat']);

%% ------------------ MAP WITH DT FOR THIS EVENT ------------------

figure(31), clf, hold on
mkfig_CascMAP
set(gcf,'position',[200 300 630 630*plot_size_ratio( lonlims,latlims )])


nnan = ~isnan(nanmean(all_dT,2));
[ sta_terms,evt_terms ] = lsq_sta_evt( all_dT,0.01 );
% hp = scatter(slons(nnan),slats(nnan),140,nanmean(all_dT(nnan,:),2),'filled','MarkerEdgeColor','k');
% hp = scatter(slons(nnan),slats(nnan),10*sum(~isnan(all_dT(nnan,:)),2),nanmean(all_dT(nnan,:),2),'filled','MarkerEdgeColor','k');
hp = scatter(slons(nnan),slats(nnan),scale*sqrt(sum(~isnan(all_dT(nnan,:)),2)),sta_terms(nnan),'filled','MarkerEdgeColor','k');

colormap(cmap), caxis(tdflim)



% % plot the raw differential travel times
% dT = [eqar.dT]';
% seaz  = [eqar.seaz]';
% slats = [eqar.slat]';
% slons = [eqar.slon]';
% 
% hold on
% for ia = 1:length(tstar)
%     line(slons(ia)+[0 sind(seaz(ia))]*scale,...
%          slats(ia)+[0 cosd(seaz(ia))]*scale,...
%          'LineWidth',1.9,...
%          'color',colour_get(tstar(ia),tstlim(2),tstlim(1)));
% end
% hold off

%% colour bar
cbar_custom(gca, 'location',[-131.4 -131.0 39.2 42.5],'tickside','right',...
    'lims',tdflim,'tickvals',[tdflim(1):ip*0.5:tdflim(2)],'cmap',cmap,...
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
title(sprintf('$\\delta T$ for $%s$-waves (%s component)',phase,component),...
      'FontSize',18,'FontWeight','bold','Interpreter','latex')
set(gca,'FontSize',14)
box on

% save
if ifsave
save2pdf(31,sprintf('dT_%s_%s_sta_average',phase,component),'figs');
end

end
