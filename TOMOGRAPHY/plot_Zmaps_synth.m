% Paper map and cross section plots in one place
clear all
HQmin = 0.2; % minimum HQ to plot
weq = 0;

def_z = 3; % default horizontal layer to show
Xmax = 15; % distance from section to plot events, in km
topo_exaggerate = 3.5;


%% Starts
wdir = '~/Documents/MATLAB/PNG_tomog/SANIS/myinv2/';
rdir = '~/Work/Papers/PNG_tomog_anisotropic/Codes_and_SVGs/tomog_slices/';
cd(wdir);


%% Bounds
min_lon = 148.5; %min(mod2.lon(ind_z))+2.5;
max_lon = 152; %max(mod2.lon(ind_z))-1.5;
min_lat = -11; %min(mod2.lat(ind_z))+2.5;
max_lat = -8; %max(mod2.lat(ind_z))-3;

%% LOAD DATA
load model
load par
grd = conv2plotable( model,par );
coast=load('~/Documents/MATLAB/PNG_tomog/papuacoast.mat');
topography = load('~/Documents/MATLAB/PNG_tomog/etopo1_bedrock.xyz');
moho = load('~/Documents/MATLAB/PNG_tomog/moho_jingle_init.mat'); %loads moho struct
load('~/Documents/MATLAB/PNG_tomog/SANIS/aniscmap.mat');

% % Get local seismicity data
% CDP = struct([]);
% dbnam = 'cdp_repick'; dbdir = '/Volumes/Zach/Documents/MATLAB/PNG_localeq/cdp_dbloc/';
% db = dbopen([dbdir dbnam],'r');
% dbor = dblookup_table(db,'origin');
% dborerr = dblookup_table(db,'origerr');
% dbjoer = dbjoin(dbor,dborerr);
% dbnrecs(dbjoer)
% dbjoers= dbsubset(dbjoer,'sdepth < 5 && ndef > 10 && smajax < 1');
% dbnrecs(dbjoers)
% dbjoers= dbsubset(dbjoers,sprintf('lat <= %f && lat >= %f && lon <= %f && lon >= %f',max_lat,min_lat,max_lon, min_lon));
% dbnrecs(dbjoers)
% [CDP(1).orid,CDP(1).elat,CDP(1).elon,CDP(1).edep,CDP(1).evtime] = dbgetv(dbjoers,'orid','lat','lon','depth','time');
% [CDP.sxx,CDP.syy,CDP.szz,CDP.smajax,CDP.sdepth] = dbgetv(dbjoers,'sxx','syy','szz','smajax','sdepth');
% dbclose(db);
% % Get EHB catalog data
% EHB = struct([]);
% dbnam = 'newehbdb'; dbdir = '~/data/catalogues/EHB/';
% db = dbopen([dbdir dbnam],'r');
% dbor = dblookup_table(db,'origin');
% dborerr = dblookup_table(db,'origerr');
% dbjoer = dbjoin(dbor,dborerr);
% dbjoers= dbsubset(dbjoer,'smajax < 20 && ndp > 0 && sdepth < 20');
% dbjoers= dbsubset(dbjoers,sprintf('lat <= %f && lat >= %f && lon <= %f && lon >= %f',max_lat,min_lat,max_lon, min_lon));
% [EHB(1).orid,EHB(1).elat,EHB(1).elon,EHB(1).edep,EHB(1).evtime] = dbgetv(dbjoers,'orid','lat','lon','depth','time');
% dbclose(db);
% 
% elat = [EHB.elat; CDP.elat];
% elon = [EHB.elon; CDP.elon];
% edep = [EHB.edep; CDP.edep];

%% GEOMETRY 
zz = unique(grd.z);
yy = [ min_lat:0.1:max_lat ]';
xx = min_lon:0.1:max_lon;

wpos=[1 181 600 530]; % ZE old window position
vanstr = {'V','An'};

%% PLOT VERTICAL SLICES
sxys(:,:,1) = [149.75 -10.9; 149.75 -8.2]; 
sxys(:,:,2) = [150.7 -10.9; 150.7 -8.2];
sxys(:,:,3) = [148.7 -9.6; 151.7 -9.6];
% sxys(:,:,4) = [150.2 -10.55; 150.2 -8.35];

for van = 0:1
	%% Velocity or Anisotropy parms
    if van==0
    val = grd.dv(:,:,def_z); %#ok<NASGU> % velocity
    cmap = flipud(jet); cbounds = [-6 6]; %#ok<NASGU>
    elseif van==1
    val = grd.va(:,:,def_z); % anisotropy
    cmap = flipud(cc); cbounds = [-4 4];
    end

    %% plot horiz slice
    figure(14), clf, set(gcf,'position',wpos), 

    hold on
    % mask poor HQ
    % val(hq<HQmin) = NaN
    contourf(grd.ln(:,:,def_z),grd.lt(:,:,def_z),val,160);
    % contour HQ
    [~,h] = contour(grd.ln(:,:,def_z),grd.lt(:,:,def_z),grd.hq(:,:,def_z),[HQmin:0.1:0.9]);
    set(h,'LineColor','k','LineStyle','--');
    
    shading flat
    daspect([ 1 cosd(abs((min_lat+max_lat)/2)) 1])
    colormap(cmap); caxis(cbounds); colorbar;
    text(min_lon + 0.1, min_lat+0.14,[num2str(zz(def_z)),' km'],'FontSize',18,'FontWeight','bold')
    geoshow(coast.ncst(:,2), coast.ncst(:,1), 'Color', 'black','linewidth',1.5)

    set(gca,'Layer','Top','FontSize',14)
    xt=[min_lon+0.5:0.5:max_lon]; yt= [min_lat:0.5:max_lat];
    set(gca,'XTick',xt,'YTick',yt) 
    set(gcf,'PaperPositionMode','auto');  
    xlim([min_lon max_lon])
    ylim([min_lat max_lat])
    box on

%     if weq
%     % find local seismicity
%         ind70 = find(edep > 80); % to plot intermediate depths EQs
%         he = plot(elon,elat,'o');
%         he2 = plot(elon(ind70),elat(ind70),'o');
%         set(he,'MarkerEdgeColor','k','MarkerFaceColor',[0.9 0.9 0.9],'MarkerSize',8,'LineWidth',1)
%         set(he2,'MarkerEdgeColor','k','MarkerFaceColor','m','MarkerSize',8,'LineWidth',1)
%     end

for ixy = 1:size(sxys,3)
    sxy = sxys(:,:,ixy);
    istr = char(64+[-1,0]+2*(ixy));
    
    figure(14); hold on
    h1 = line(sxy(:,1),sxy(:,2)); 
    set(h1,'LineWidth',2.5,'LineStyle','--','color','k')
    if diff(sxy,2) > 0 % N-S section
    text(sxy(1,1)-0.1,sxy(1,2)-0.02,istr(1),'HorizontalAlignment','center','FontWeight','bold','FontSize',15)
    text(sxy(2,1)-0.1,sxy(2,2)+0.02,istr(2),'HorizontalAlignment','center','FontWeight','bold','FontSize',15)
    else
    text(sxy(1,1)-0.05,sxy(1,2),istr(1),'HorizontalAlignment','right','FontWeight','bold','FontSize',15)
    text(sxy(2,1)+0.05,sxy(2,2),istr(2),'HorizontalAlignment','left','FontWeight','bold','FontSize',15)
    end


% plot section
% conditions are to get nice looking plot contextualised geographically
figure(15), clf, set(gcf,'position',wpos + [700 0 0 0]), 
hold on

% pick section

[dV,anis,hqz,slt,sln,ss,Dxy] = get_depth_section(sxy(1,:),sxy(2,:),grd);
if van==0, val=dV; elseif van==1; val=anis; end
Xxy = distance_km(sxy(1,2),sxy(1,1),sxy(2,2),sxy(2,1));
ss = linspace(0,Xxy,length(ss));
ss = ss-mean(ss);

% % mask out low hq
% dVz(hqz<HQmin) = NaN;
% hqz(zz<40,:)=NaN;


%Nicer version: plot against distance - equal axes
contourf(ss,zz,val,160);
set(gca,'YDir','reverse','FontSize',14);
axis equal

shading flat
ylabel('Depth (km)','FontSize',16)
xlabel('Distance along profile (km)','FontSize',16)
colormap(cmap); caxis(cbounds)

% contour HQ
[cs,h] = contour(ss,zz,hqz,[HQmin:0.1:0.9]);
set(h,'LineColor','k','LineStyle','--');

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

% plot on moho
% mohh = interp2(moho.xi,moho.yi,moho.moho_surf,sxx,syy); % if using moho_pick
mohh = interp2(moho.lon,moho.lat,moho.mohodepth,sln,slt); % if using moho_jingle_*
plot(ss,mohh,'-k','LineWidth',3)

% plot crust
mhn1 = find(~isnan(mohh),1,'first'); mohh(1:mhn1-1)   = mohh(mhn1);
mhn2 = find(~isnan(mohh),1,'last');  mohh(mhn2+1:end) = mohh(mhn2);

crust.x = [ss(1)-50,ss,ss(end)+50,ss(end)+50,fliplr(ss),ss(1)-50,ss(1)-50]';
crust.y = [topo(1);topo;topo(end);mohh(end);flipud(mohh);mohh(1);topo(1)];
patch(crust.x,crust.y,[0.85 0.85 0.5],'LineStyle',':','LineWidth',0.5);

% % re-lay topog & moho
% plot(ss,mohh,'-k','LineWidth',2)
% plot(ss,topo,'k','LineWidth',2); 


% plot box and axes
set(gca,'Layer','Top','FontSize',14,...
    'YTick',[-shifter,5,50:50:250]','YTickLabel',num2str([0,5,50:50:250]'))
ylim([-shifter - 3*topo_exaggerate,zz(end-1)])
xlim([min(ss), max(ss)])
box on

% plot section ends
text(min(ss),-30,istr(1),'HorizontalAlignment','center','FontWeight','bold','FontSize',16)
text(max(ss),-30,istr(2),'HorizontalAlignment','center','FontWeight','bold','FontSize',16)

if weq
% plot local seismicity
% calc. distance of earthquakes to line, which are in bounds
Dev = dist2line(sxy(1,:),sxy(2,:),[elon,elat])*Xxy/norm(diff(sxy));
indx = find(abs(Dev)<Xmax);
% calc. projection of each earthquakes along line
Sev = ([elon,elat]-ones(length(elon),1)*sxy(1,:))*diff(sxy)'*Xxy/(norm(diff(sxy))^2);
% only include earthquakes within end-bounds of line
indx1 = find(Sev<=Xxy & Sev>=0);
indx = intersect(indx,indx1);

he = plot(Sev(indx)-(Xxy/2),edep(indx),'o');
set(he,'MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',8,'MarkerFaceColor','w')
end

set(gcf, 'Color', [1 1 1])
save2pdf(15,strcat(vanstr{van+1},'_Zmap_',num2str(ixy)),rdir);
% export_fig '/Users/Zach/Documents/MATLAB/PNG_tomog/figs/Zslice_weq' -eps -transparent -r275
% movefile('~/Documents/MATLAB/PNG_tomog/figs/Zslice_weq.eps',['~/Documents/MATLAB/PNG_tomog/figs/Zslice_weq' num2str(ixy) '.eps']);

end

save2pdf(14,strcat(vanstr{van+1},'_Zmap_horiz'),rdir);

end


