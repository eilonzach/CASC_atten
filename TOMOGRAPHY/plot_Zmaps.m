function plot_Zmaps(plot_model,par,data,ifsave)

HQmin = 0.3; % minimum HQ to plot

def_z = 3; % default horizontal layer to show

topo_exaggerate = 3.5;

%% PLOT VERTICAL SLICES
sxys(:,:,1) = [-132 47; -118 46.5]; 
sxys(:,:,2) = [-132 44; -118 44]; 
sxys(:,:,3) = [-132 42; -118 41]; 
sxys(:,:,4) = [-121.5 39.5; -121.5 49];

%% Bounds
min_lon = -132; %min(mod2.lon(ind_z))+2.5;
max_lon = -117.5; %max(mod2.lon(ind_z))-1.5;
min_lat = 39; %min(mod2.lat(ind_z))+2.5;
max_lat = 50; %max(mod2.lat(ind_z))-3;

%% Velocity or Q parms
if par.t_ts==1
cmp = flipud(jet);
cbounds = [-3 3];
elseif par.t_ts==2
cmp = parula; 
cbounds = [-600 200];
end

%% Starts
wdir = '~/Documents/MATLAB/CASC_atten/TOMOGRAPHY/';
figdir = '~/Documents/MATLAB/CASC_atten/TOMOGRAPHY/figs/';
grd = plot_model;
wpos=[1 500 1000 530]; % window position
typestr = {'dV','dq'};

%% LOAD DATA
cd('~/Documents/MATLAB/CASC_atten/mapdata/');
coast=load('m_casccoast.mat');
[topoX,topoY,topoZ] = grdread2('cascmrg2.grd'); % load topo grid
[topoX,topoY] = meshgrid(double(topoX),double(topoY)); topoZ = double(topoZ);
jdf = dlmread('ridge_xy'); % load ridge
fzs = dlmread('transforms_xy'); % load transforms & fracture zones
[Vnam,Vlon,Vlat,~,~] = textread('volcanoes_complete','%s %f %f %s %s','headerlines',1);
cd(wdir);

%% GEOMETRY 
zz = unique(grd.z);
yy = [ min_lat:0.1:max_lat ]';
xx = min_lon:0.1:max_lon;


%% plot horiz slice
figure(14), clf, set(gcf,'position',[1 500 590 750]), hold on

% mask poor HQ
valz(grd.hq(:,:,def_z)<HQmin) = NaN;
contourf(grd.ln(:,:,def_z),grd.lt(:,:,def_z),100*grd.val(:,:,def_z),160);
% contour HQ
[~,h] = contour(grd.ln(:,:,def_z),grd.lt(:,:,def_z),grd.hq(:,:,def_z),[HQmin:0.1:0.9]);
set(h,'LineColor','k','LineStyle','--','LineWidth',0.7);

shading flat
daspect([ 1 cosd(abs((min_lat+max_lat)/2)) 1])
colormap(cmp); caxis(cbounds); colorbar;
text(min_lon + 0.5, min_lat+0.5,[num2str(zz(def_z)),' km'],'FontSize',22,'FontWeight','bold')
geoshow(coast.ncst(:,2), coast.ncst(:,1), 'Color', 'black','linewidth',2)
geoshow(jdf(:,2), jdf(:,1), 'Color', 'black','linewidth',1)
geoshow(fzs(:,2), fzs(:,1),'Linestyle','--', 'Color', 'black','linewidth',1)
    

set(gca,'Layer','Top','FontSize',16)
xt=[min_lon:2:max_lon]; yt= [min_lat:2:max_lat];
set(gca,'XTick',xt,'YTick',yt) 
set(gcf,'PaperPositionMode','auto');  
xlim([min_lon max_lon])
ylim([min_lat max_lat])
box on


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
    figure(14+ixy), clf, set(gcf,'position',wpos + [700 0 0 0]), 
    hold on

    %% get model along the section
    [valz,hqz,slt,sln,ss,Dxy,smb,Xxy] = get_depth_section(sxy(1,:),sxy(2,:),grd);

    ss = linspace(0,Xxy,length(ss));
    ss = ss-mean(ss);

    % % mask out low hq
    valz(hqz<HQmin) = NaN;
    % hqz(zz<40,:)=NaN;

    contourf(ss,zz,100*valz,160);
    set(gca,'YDir','reverse','FontSize',16);
    axis equal

    shading flat
    ylabel('Depth (km)','FontSize',20,'FontWeight','bold')
    xlabel('Distance along profile (km)','FontSize',20,'FontWeight','bold')
    colormap(cmp); caxis(cbounds)

    % contour HQ
    [~,h] = contour(ss,zz,hqz,[0.6:0.1:0.9]);
    set(h,'LineColor','k','LineStyle','--','LineWidth',0.7);

    % % contour semb
    % [cs,h] = contour(ss,zz,smb,0.7);
    % set(h,'LineColor','b','LineStyle','--','LineWidth',2);
    

    % CHANGE SCALE AT 5 KM DEPTH
    shifter = 5*(topo_exaggerate - 1); % from simple maths

    %% plot on topography
    topo = interp2(topoX,topoY,topoZ,sln,slt);
    topo = -topo./(1000./topo_exaggerate);
        % CHANGE SCALE AT 5 KM DEPTH
    topo = topo-shifter;
    ht = plot(ss,topo,'k','LineWidth',3);

    % plot water
    water.x = [ss,fliplr(ss)]';
    % CHANGE SCALE AT 5 KM DEPTH (otherwise would be water.y = [zeros(size(topo)); flipud(topo)]; water.y(water.y<0)=0;
    water.y = [-shifter*ones(size(topo)); flipud(topo)]; water.y(water.y<-shifter)=-shifter;
    hw = patch(water.x,water.y,'b');
    set(hw,'EdgeAlpha',0,'FaceAlpha',0.5)
    
    %% plot on slab
    [slbs,slbz] = slab_on_section(sxy);
    hsl = plot(slbs-(Xxy/2),-slbz,'--k','LineWidth',3);

    
    %% plot volcanoes within 50 km
    % calc. distance of stations to line, which are in bounds
    Dev = dist2line(sxy(1,:),sxy(2,:),[Vlon,Vlat])*Xxy/norm(diff(sxy));
    indx = find(abs(Dev)<50);
    % calc. projection of each stations along line
    Sev = ([Vlon,Vlat]-ones(length(Vlon),1)*sxy(1,:))*diff(sxy)'*Xxy/(norm(diff(sxy))^2);
    % only include stations within end-bounds of line
    indx1 = find(Sev<=Xxy & Sev>=0);
    indx = intersect(indx,indx1);
    hs = plot(Sev(indx)-(Xxy/2),interp1(ss,topo,Sev(indx)-(Xxy/2)),'^');
    set(hs,'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','r')
    
    %% PLOT STATIONS within 50 km
    % calc. distance of stations to line, which are in bounds
    Dev = dist2line(sxy(1,:),sxy(2,:),[data.stn.lon,data.stn.lat])*Xxy/norm(diff(sxy));
    indx = find(abs(Dev)<50);
    % calc. projection of each stations along line
    Sev = ([data.stn.lon,data.stn.lat]-ones(length(data.stn.lon),1)*sxy(1,:))*diff(sxy)'*Xxy/(norm(diff(sxy))^2);
    % only include stations within end-bounds of line
    indx1 = find(Sev<=Xxy & Sev>=0);
    indx = intersect(indx,indx1);
    hs = plot(Sev(indx)-(Xxy/2),-data.stn.elv(indx)*topo_exaggerate-shifter,'v');
    set(hs,'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','c')
    uistack(hs,'top');
    uistack(hw,'bottom');
    uistack(ht,'bottom')
    
    
    
    %% plot box and axes
    set(gca,'Layer','Top','FontSize',16,...
        'YTick',[-shifter,5,50:50:400]','YTickLabel',num2str([0,5,50:50:400]'))
    ylim([-shifter - 5*topo_exaggerate,zz(end-1)])
    xlim([min(ss), max(ss)])
    box on

    % plot section ends
    text(min(ss),-33,istr(1),'HorizontalAlignment','center','FontWeight','bold','FontSize',24)
    text(max(ss),-33,istr(2),'HorizontalAlignment','center','FontWeight','bold','FontSize',24)

    %% SAVE
    set(gcf, 'Color', [1 1 1])
    if ifsave
        save2pdf(14+ixy,strcat(typestr{par.t_ts},'_Zmap',num2str(ixy)),figdir);
    end

end % loop on sections

%% Save hmap
    if ifsave
        save2pdf(14,strcat(typestr{par.t_ts},'_Zmap_hplot'),figdir);
    end

end

function [slbs,slbz] = slab_on_section(sxy)
    Xxy = distance_km(sxy(1,2),sxy(1,1),sxy(2,2),sxy(2,1));
    [ slab_contrs ] = McCrory_slab_depth;
    Nc = length(slab_contrs);
    slbs = nan(Nc,1);
    slbz = nan(Nc,1);
    for ic = 1:Nc
        dst = dist2line(sxy(1,:),sxy(2,:),[slab_contrs(ic).lon,slab_contrs(ic).lat]);
        if ~any(dst<0), continue; end
        if ~any(dst>0), continue; end
        ind = [find(dst==min(dst(dst>0))),find(dst==max(dst(dst<0)))];        
        ss = ([slab_contrs(ic).lon(ind),slab_contrs(ic).lat(ind)]-ones(2,1)*sxy(1,:))*diff(sxy)'*Xxy/(norm(diff(sxy))^2);
        slbs(ic) = mean(ss);
        slbz(ic) = slab_contrs(ic).depth;
    end
    
    % fill in any nans w/ linear interp
    nnan = ~isnan(slbz);
    slbs = interp1(slbz(nnan),slbs(nnan),[slab_contrs.depth]);
    slbz = interp1(slbz(nnan),slbz(nnan),[slab_contrs.depth]);
end

