function plot_basic(plot_model,par,opt,saveopt)
% opt is an option describing the input model
%  opt == 1  means real result
%  opt == 2  means synth in
%  opt == 3  means synth out
%  opt == 4  means plot bootstrap uncertainties
%
% saveopt is an option to save (1) or not (0)
if nargin < 3
    opt = 1;
end
if nargin < 4
    saveopt = 0;
end

clims = 5*[-1 1];

if par.t_ts == 1
    cmp = flipud(jet);
    valstr = 'V';
elseif par.t_ts == 2;
    cmp = parula;
    valstr = 'q';
end

min_lon = -131; %min(mod2.lon(ind_z))+2.5;
max_lon = -120; %max(mod2.lon(ind_z))-1.5;
min_lat = 39; %min(mod2.lat(ind_z))+2.5;
max_lat = 51; %max(mod2.lat(ind_z))-3;

%% load some data
mapdata = '/Users/zeilon/Documents/MATLAB/CASC_atten/mapdata/';
coast = load([mapdata,'m_casccoast.mat']); % load coastline
jdf = dlmread([mapdata,'ridge_xy']); % load ridge
fzs = dlmread([mapdata,'transforms_xy']); % load transforms & fracture zones


%% ===========================  PLOT VELOCITY  ===========================
%% ===========================  PLOT VELOCITY  ===========================
figure(57 + 10*opt); clf

set(gcf,'Position',[0,280-20*opt,900,700])
for iz = 2:par.nz-1
    subplot(3,4,iz-1)
	hold on
    
    yy = [ min_lat:0.1:max_lat ]';
    xx = min_lon:0.1:max_lon;

    if opt~=4
        VAL = griddata(plot_model.ln(:,:,iz),plot_model.lt(:,:,iz),100*plot_model.val(:,:,iz),xx,yy);        
        cbounds = clims;
    elseif opt==4
        VAL = griddata(plot_model.ln(:,:,iz),plot_model.lt(:,:,iz),100*plot_model.sv(:,:,iz),xx,yy);
        cbounds = [0 0.1];
        cmp = colormap(cmap_makecustom([0.8 0.8 0.1],[0.1 0.8 0.8],0));
    end
%     smb = griddata(plot_model.ln(:,:,iz),plot_model.lt(:,:,iz),plot_model.semb_v(:,:,iz),xx,yy);
    hq = griddata(plot_model.ln(:,:,iz),plot_model.lt(:,:,iz),plot_model.hq(:,:,iz),xx,yy);
    
    contourf(xx,yy,VAL,80);
%     contour(xx,yy,smb,[0.7:0.1:1],'--r','Linewidth',1.5)
%     contour(xx,yy,smb.*hq,[0.5:0.1:1],'--b','Linewidth',1.5)
    
    shading flat
    axis( [min_lon max_lon min_lat max_lat] );
    daspect([ 1 cosd(abs((min_lat+max_lat)/2)) 1])
    colormap(cmp)
    caxis(cbounds)
    
    geoshow(coast.ncst(:,2), coast.ncst(:,1), 'Color', 'black','linewidth',2)
    geoshow(jdf(:,2), jdf(:,1), 'Color', 'black','linewidth',1)
    geoshow(fzs(:,2), fzs(:,1),'Linestyle','--', 'Color', 'black','linewidth',1)
    
    xlabel('Longitude'); 
    ylabel('Latitude');
    title(sprintf('Depth slice %.0f km',par.zz(iz)),'FontSize',14)
end

%% scale


subplot(3,4,iz)
set(gca,'Visible','off')
if opt~=4
    str = ['$\delta ',valstr,' \,\, \%$']; 
    colormap(cmp); 
    caxis(clims)
elseif opt==4
    str = '%\sigma_{boot}$'; 
    colormap(cmap_makecustom([0.8 0.8 0.1],[0.1 0.8 0.8],0)); 
    caxis([0 0.1]);
end
hc = colorbar('EastOutside');
set(get(hc,'YLabel'),'String',str,'FontSize',20,'FontWeight','bold','interpreter','latex')
set(hc,'FontSize',12,'FontWeight','bold','Position',[0.8 0.12 0.02 0.2])


%% save
if saveopt
    suff = {'all','all_synin','all_synout','booterrs'};
    pref = ['d',valstr];
    if par.wtdata == 0, wtstr = '_nowt'; else wtstr = '';  end
   
    fprintf('Saving figure %s_%s%s... \n',pref,suff{opt},wtstr);
    
    ostr = sprintf('figs/%s_%s%s.eps',pref,suff{opt},wtstr);
    print(ostr,'-depsc','-r600');

end

end