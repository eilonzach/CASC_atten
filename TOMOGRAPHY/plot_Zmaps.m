function plot_Zmaps(plot_model,par)

HQmin = 0.5; % minimum HQ to plot

def_z = 4; % default horizontal layer to show
Xmax = 15; % distance from section to plot events, in km
topo_exaggerate = 3.5;

%% PLOT VERTICAL SLICES
sxys(:,:,1) = [149.8 -10.9; 149.8 -8.2]; 
sxys(:,:,2) = [150.4 -10.9; 150.4 -8.2];
sxys(:,:,3) = [151 -10.9; 151 -8.2];
% sxys(:,:,3) = [151 -8.7; 151.7 -8.7];
% sxys(:,:,4) = [150.2 -10.55; 150.2 -8.35];

grd = plot_model;

%% Velocity or Q parms
if van==0
val = grd.dv(:,:,def_z); %#ok<NASGU> % velocity
cmap = flipud(jet); cbounds = cbounds_v; %#ok<NASGU>
cbounds = [-6 6];
elseif van==1
val = grd.va(:,:,def_z); % anisotropy
cmap = flipud(cc); cbounds = cbounds_a;
cbounds = [-6 6];
end

%% Starts
wdir = '~/Documents/MATLAB/PNG_tomog/SANIS/myinv2/';
figdir = '~/Work/Papers/PNG_tomog_anisotropic/Codes_and_SVGs/tomog_slices/';
cd(wdir);


%% Bounds
min_lon = -132; %min(mod2.lon(ind_z))+2.5;
max_lon = -120; %max(mod2.lon(ind_z))-1.5;
min_lat = 39; %min(mod2.lat(ind_z))+2.5;
max_lat = 50; %max(mod2.lat(ind_z))-3;

%% LOAD DATA
coast=load('~/Documents/MATLAB/CASC_atten/mapdata/m_casccoast.mat');
[grdX,grdY,grdZ] = grdread2('~/Documents/MATLAB/CASC_atten/mapdata/cascmrg2.grd'); % load topo grid
[grdX,grdY] = meshgrid(double(grdX),double(grdY)); grdZ = double(grdZ);


%% GEOMETRY 
zz = unique(grd.z);
yy = [ min_lat:0.1:max_lat ]';
xx = min_lon:0.1:max_lon;

wpos=[1 500 600 530]; % ZE old window position
vanstr = {'V','An'};



%% plot horiz slice
figure(14), clf, set(gcf,'position',wpos), 

hold on
% mask poor HQ
val(grd.hq(:,:,def_z)<HQmin) = NaN;
contourf(grd.ln(:,:,def_z),grd.lt(:,:,def_z),val,160);
% contour HQ
[~,h] = contour(grd.ln(:,:,def_z),grd.lt(:,:,def_z),grd.hq(:,:,def_z),[HQmin:0.1:0.9]);
set(h,'LineColor','k','LineStyle','--');

shading flat
daspect([ 1 cosd(abs((min_lat+max_lat)/2)) 1])
colormap(cmap); caxis(cbounds); colorbar;
text(min_lon + 0.1, min_lat+0.14,[num2str(zz(def_z)),' km'],'FontSize',22,'FontWeight','bold')
geoshow(coast.ncst(:,2), coast.ncst(:,1), 'Color', 'black','linewidth',1.5)

set(gca,'Layer','Top','FontSize',16)
xt=[min_lon+0.5:0.5:max_lon]; yt= [min_lat:0.5:max_lat];
set(gca,'XTick',xt,'YTick',yt) 
set(gcf,'PaperPositionMode','auto');  
xlim([min_lon max_lon])
ylim([min_lat max_lat])
box on

if weq
% find local seismicity
    ind70 = find(edep > 80); % to plot intermediate depths EQs
    he = plot(elon,elat,'o');
    he2 = plot(elon(ind70),elat(ind70),'o');
    set(he,'MarkerEdgeColor','k','MarkerFaceColor',[0.9 0.9 0.9],'MarkerSize',8,'LineWidth',1)
    set(he2,'MarkerEdgeColor','k','MarkerFaceColor','m','MarkerSize',8,'LineWidth',1)
end

for ixy = 1:size(sxys,3)
    sxy = sxys(:,:,ixy);
    istr = char(64+[-1,0]+2*(ixy));
    
    figure(14); hold on
    h1 = line(sxy(:,1),sxy(:,2)); 
    set(h1,'LineWidth',2.5,'LineStyle','--','color','k')
    if diff(sxy,2) > 0 % N-S section
    text(sxy(1,1)-0.1,sxy(1,2)-0.02,istr(1),'HorizontalAlignment','center','FontWeight','bold','FontSize',22)
    text(sxy(2,1)-0.1,sxy(2,2)+0.02,istr(2),'HorizontalAlignment','center','FontWeight','bold','FontSize',22)
    else
    text(sxy(1,1)-0.05,sxy(1,2),istr(1),'HorizontalAlignment','right','FontWeight','bold','FontSize',22)
    text(sxy(2,1)+0.05,sxy(2,2),istr(2),'HorizontalAlignment','left','FontWeight','bold','FontSize',22)
    end


% plot section
% conditions are to get nice looking plot contextualised geographically
figure(15), clf, set(gcf,'position',wpos + [700 0 0 0]), 
hold on

% pick section

[dV,anis,hqz,slt,sln,ss,Dxy,smbv,smba,Xxy] = get_depth_section(sxy(1,:),sxy(2,:),grd);
if van==0, val=dV; elseif van==1; val=anis; end
if van==0, smb=smbv; elseif van==1; smb=smba; end
ss = linspace(0,Xxy,length(ss));
ss = ss-mean(ss);

% % mask out low hq
val(hqz<HQmin) = NaN;
% hqz(zz<40,:)=NaN;


%Nicer version: plot against distance - equal axes
contourf(ss,zz,val,160);
set(gca,'YDir','reverse','FontSize',16);
axis equal

shading flat
ylabel('Depth (km)','FontSize',20,'FontWeight','bold')
xlabel('Distance along profile (km)','FontSize',20,'FontWeight','bold')
colormap(cmap); caxis(cbounds)

% contour HQ
[cs,h] = contour(ss,zz,hqz,[0.6:0.1:0.9]);
set(h,'LineColor','k','LineStyle','--');

% % contour semb
% [cs,h] = contour(ss,zz,smb,0.7);
% set(h,'LineColor','b','LineStyle','--','LineWidth',2);

% CHANGE SCALE AT 5 KM DEPTH
shifter = 5*(topo_exaggerate - 1); % from simple maths

% plot on topography
topo = griddata(topography(:,1),topography(:,2),topography(:,3),sln,slt);
topo = -topo./(1000./topo_exaggerate);
    % CHANGE SCALE AT 5 KM DEPTH
topo = topo-shifter;
plot(ss,topo,'k','LineWidth',3)
 
% plot water
water.x = [ss,fliplr(ss)]';
% CHANGE SCALE AT 5 KM DEPTH (otherwise would be water.y = [zeros(size(topo)); flipud(topo)]; water.y(water.y<0)=0;
water.y = [-shifter*ones(size(topo)); flipud(topo)]; water.y(water.y<-shifter)=-shifter;
hw = patch(water.x,water.y,'b');
set(hw,'EdgeAlpha',0,'FaceAlpha',0.5)


% % re-lay topog & moho
% plot(ss,mohh,'-k','LineWidth',2)
% plot(ss,topo,'k','LineWidth',2); 


% plot box and axes
set(gca,'Layer','Top','FontSize',16,...
    'YTick',[-shifter,5,50:50:250]','YTickLabel',num2str([0,5,50:50:250]'))
ylim([-shifter - 3*topo_exaggerate,zz(end-1)])
ylim([-shifter - 3*topo_exaggerate,280])
xlim([min(ss), max(ss)])
box on

% plot section ends
text(min(ss),-33,istr(1),'HorizontalAlignment','center','FontWeight','bold','FontSize',24)
text(max(ss),-33,istr(2),'HorizontalAlignment','center','FontWeight','bold','FontSize',24)

return

set(gcf, 'Color', [1 1 1])
save2pdf(15,strcat(vanstr{van+1},'_Zmap_',num2str(ixy)),rdir);

save2pdf_old(14,[vanstr{van+1},'_Zmap_horiz',num2str(zz(def_z))],rdir);

end % loop on sections


end



