%% Make a nice map of the Cascadia region, as a basemap for further plotting

addpath('/Users/zeilon/Documents/MATLAB/lib/grdread2')
addpath('/Users/zeilon/Documents/MATLAB/lib/haxby_colourmap')

if ~exist('latlim','var'), latlim=[39 52];      end
if ~exist('lonlim','var'), lonlim=[-133 -120];  end

cd('/Users/zeilon/Documents/MATLAB/CascAtten/mapdata');

coast = load('m_casccoast.mat'); % load coastline
[grdX,grdY,grdZ] = grdread2('cascmrg2.grd'); % load topo grid
[grdX,grdY] = meshgrid(double(grdX),double(grdY)); grdZ = double(grdZ);

res2 = dlmread('2D_bathym_section_2.txt','\t',1,0); % load section
res3 = dlmread('2D_bathym_section_3.txt','\t',1,0); % load section
% res = dlmread('2D_bathym_section.txt','\t',1,0); % load section
res = [[res2(:,1:3);res3(:,1:3)],[res2(:,4);res2(end,4)+res3(:,4)]];
    section_lola = res(:,1:2);
    section_z = res(:,3);
    section_x = res(:,4);
    
    
jdf = dlmread('ridge_xy'); % load ridge

fzs = dlmread('transforms_xy'); % load transforms & fracture zones

[Vnam,Vlon,Vlat,~,~] = textread('volcanoes','%s %f %f %s %s','headerlines',1);


clf, set(gcf,'position',[200 200 600 800])

%% Un-comment to make a new coastline object
% m_proj('mercator','longitudes',lonlim,'latitudes',latlim);
% m_gshhs_h('patch',[.7 .7 .7],'edgecolor','k');
% m_gshhs_h('save','m_casccoast');
% m_usercoast('m_casccoast','patch',[.7 .7 .7],'edgecolor','k');
% m_coord('geographic');
% m_grid;

hold on


%% make basemap image
image(grdX(:),grdY(:),colour_get(grdZ,3000,-5000,haxby)); % where numbers are maxZ, minZ, cmap

plot(coast.ncst(:,1),coast.ncst(:,2),'k','LineWidth',1)

scatter( Vlon, Vlat,65,'k','^','filled') %Volcanoes

plot(jdf(:,1),jdf(:,2),'k','Linewidth',2) % ridge
plot(fzs(:,1),fzs(:,2),'--k','Linewidth',1) % transforms

% plot section
plot(section_lola(:,1),section_lola(:,2),'k','LineWidth',2)
plot(linterp(section_x,section_lola(:,1),[-100:100:700]'),...
     linterp(section_x,section_lola(:,2),[-100:100:700]'),...
     '.k','MarkerSize',25)
text(linterp(section_x,section_lola(:,1),[-100:100:700]'),...
     linterp(section_x,section_lola(:,2),[-100:100:700]')+0.12,...
     num2str([-100:100:700]'),'FontWeight','bold')


xlim(lonlim)
ylim(latlim)
grid on
box on

cd ..