%% Make a nice map of the Cascadia region, as a basemap for further plotting

addpath('/Users/zeilon/Documents/MATLAB/lib/grdread2')
addpath('/Users/zeilon/Documents/MATLAB/lib/haxby_colourmap')
addpath('/Users/zeilon/Documents/MATLAB/CASC_atten/matguts')

if ~exist('latlims','var'), latlims=[39 51];      end
if ~exist('lonlims','var'), lonlims=[-132 -110];  end

wd = pwd;
cd('/Users/zeilon/Documents/MATLAB/CASC_atten/mapdata');

coast = load('m_casccoast.mat'); % load coastline
[grdX,grdY,grdZ] = grdread2('cascmrg2.grd'); % load topo grid
[grdX,grdY] = meshgrid(double(grdX),double(grdY)); grdZ = double(grdZ);

res2 = dlmread('2D_bathym_section_leg1.txt','\t',1,0); % load section
res3 = dlmread('2D_bathym_section_leg2.txt','\t',1,0); % load section
% res = dlmread('2D_bathym_section.txt','\t',1,0); % load section
res = [[res2(:,1:3);res3(:,1:3)],[res2(:,4);res2(end,4)+res3(:,4)]];
    section_lola = res(:,1:2);
    section_z = res(:,3);
    section_x = res(:,4);

% ofile = '2D_bathym_section_overall.txt';
% fid = fopen(ofile,'w');
% fprintf(fid,'Longitude	Latitude	Elevation, m	Distance (km)\n');
% for ii = 1:length(res);
%     fprintf(fid,'%-10.5f  %-10.6f  %-10.4f  %-11.6f\n',res(ii,:));
% end
% fclose(fid);

    
    
jdf = dlmread('ridge_xy'); % load ridge

fzs = dlmread('transforms_xy'); % load transforms & fracture zones


dx = distance_km(latlims(1),lonlims(1),latlims(1),lonlims(2));
dy = distance_km(latlims(1),lonlims(1),latlims(2),lonlims(1));

clf, set(gcf,'position',[200 200 700*[0.8*(dx/dy) 1] ])

%% Un-comment to make a new coastline object
% m_proj('mercator','longitudes',lonlims,'latitudes',latlims);
% m_gshhs_h('patch',[.7 .7 .7],'edgecolor','k');
% m_gshhs_h('save','m_casccoast');
% m_usercoast('m_casccoast','patch',[.7 .7 .7],'edgecolor','k');
% m_coord('geographic');
% m_grid;

hold on

%% make basemap image
image(grdX(:),grdY(:),colour_get(grdZ,3000,-5000,haxby)); % where numbers are maxZ, minZ, cmap

plot(coast.ncst(:,1),coast.ncst(:,2),'k','LineWidth',1)

% volcanoes
[Vnam,Vlon,Vlat,~,~] = textread('volcanoes_complete','%s %f %f %s %s','headerlines',1);
scatter( Vlon, Vlat,60,'k','^','MarkerFaceColor',[1 1 1],'LineWidth',1) %Volcanoes
[Vnam,Vlon,Vlat,~,~] = textread('volcanoes_principal','%s %f %f %s %s','headerlines',1);
scatter( Vlon, Vlat,85,'k','^','MarkerFaceColor','r','LineWidth',1) %Volcanoes


plot(jdf(:,1),jdf(:,2),'k','Linewidth',2) % ridge
plot(fzs(:,1),fzs(:,2),'--k','Linewidth',1) % transforms

%% plot section
% plot(section_lola(:,1),section_lola(:,2),'k','LineWidth',2)
% plot(linterp(section_x,section_lola(:,1),[-100:100:700]'),...
%      linterp(section_x,section_lola(:,2),[-100:100:700]'),...
%      '.k','MarkerSize',25)
% text(linterp(section_x,section_lola(:,1),[-100:100:700]'),...
%      linterp(section_x,section_lola(:,2),[-100:100:700]')+0.12,...
%      num2str([-100:100:700]'),'FontWeight','bold')

%% plot isochrons
% [chgrdX,chgrdY] = meshgrid(linspace(lonlims(1),-123,50),linspace(40,latlims(2),60)); 
% [ chgrdage,chrons,F ] = jdf_crust_age(chgrdY,chgrdX);
% % chtick = unique(chrons.age(~isnan(chrons.age))); % plot chrons
% chtick = 1:12; % plot Ma
% contour(chgrdX,chgrdY,chgrdage,chtick,'LineWidth',2,'linestyle','--')
% colormap(hot)

xlim(lonlims)
ylim(latlims)
grid on
box on

cd(wd)