function plot_delaytimes_N_E_dT(data,par,saveres)
% plot_delaytimes_N_E_dT(data,par,saveres)

figure(32), clf, hold on
mkfig_CascMAP

odir = '~/Documents/MATLAB/CASC_attem/TOMOGRAPHY/figs/';

lonlims = [-132.1 -120];
latlims = [39.1 50.9];
plotsize = 800;
dlim = [-1.6 1.6];
cmap = parula;

scale = 0.5; % length of lines

% lagtick = 1;

indx = abs(data.ray.d) < 4; 
d = data.ray.d(indx);
baz = data.ray.seaz(indx);
sd = data.ray.sd(indx);
slat = data.ray.stalat(indx);
slon = data.ray.stalon(indx);
% incs = rayp2inc(data.ray.pd(indx),4.34,115);

%% plot uncorrected times on Map
set(gcf,'position',[200 300 plotsize/plot_size_ratio(lonlims,latlims) plotsize])
set(gca,'xlim',lonlims,'ylim',latlims)


% title(sprintf('Raw time residuals (%s)',secstrs{ic}),'FontSize',18)
for ia=1:sum(indx)

% plot the data
p1la = slat(ia);
p1lo = slon(ia);
[p2la,p2lo] = reckon(p1la,p1lo,scale,baz(ia));
% plot line toward baz, coloured by tstar
plot([p1lo;p2lo],[p1la;p2la],'LineWidth',0.5/sd(ia),...
    'color',colour_get(d(ia),dlim(2),dlim(1),cmap))


end % loop on arrivals
return

hcb = colorbar('peer',gca);
set(hcb,'YTick',[0:0.5*lagtick./laglim(2):1],'YTickLabel',num2str([laglim(1):lagtick:laglim(2)]'),'FontSize',18);
set(get(hcb,'Ylabel'),'string','$\longleftarrow$  FAST ~~~~~  Tdiff (s) ~~~~~ SLOW  $\longrightarrow$',...
    'Fontsize',24,'FontWeight','bold','interpreter','Latex');
set(get(hcb,'Ylabel'),'position',get(get(hcb,'Ylabel'),'position')+[0.4 0 0])


% plot on a key of weights
wt_eg = [1, 5, 10];
sd_eg = (10*wt_eg).^-0.5;
wd_eg= 0.5./sd_eg;
lon_key = 149.05;
lat_key = -8.4;
hold on
for ii = 1:length(wt_eg)
    m_line(lon_key + [0, scale],(lat_key + 0.1*ii)*[1 1],'Linewidth',wd_eg(ii),'color','k')
    m_text(lon_key + 1.5*scale,   lat_key + 0.1*ii,num2str(wt_eg(ii)),'Fontsize',15,'FontWeight','bold')
end
% plot key with datatype
hold on
m_text(lon_key,-10.85,['\textbf{',secstrs{ic},'}'],'Fontsize',22,'FontWeight','bold','interpreter','Latex')

% save
if saveres
    ofile = sprintf('Tdiff_%s_raw',secstrs{ic});
    save2pdf(11,ofile,odir);
    copyfile([odir ofile '.pdf'],[rdir ofile '.pdf'])
end


end

