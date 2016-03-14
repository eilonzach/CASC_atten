function eqar = plot_obs1evt_COMB( eqar)

datwind = [-70 100];
filtfs = [0.05 2]; % in hz
tstlim = 1.5*[-1 1];
% refsta = 'LON';
amp_seis = 1.5; % amplification factor for the traces

t = tauptime('p',eqar(1).phase,'deg',eqar(1).gcarc); 
evtime = eqar(1).pred_arrT - t(1).time;

% only use "good" stations
gd=[]; 
for is = 1:length(eqar), 
    if ~isempty(eqar(is).dtstar_comb) && ~isnan(eqar(is).dtstar_comb), 
        gd = [gd;is]; 
    end
end
eqar = eqar(gd);
% put NaNs in for stations without specR dtstar:
for is=1:length(eqar), if isempty(eqar(is).dtstar), eqar(is).dtstar=NaN; end; end

%% get deets
phase = eqar(1).phase;
comp = eqar(1).par_dtstar_specR.comp;
window = [-eqar(1).par_dtstar_comb.wind.prex,eqar(1).par_dtstar_comb.wind.postx];
fmax = nanmean([eqar.fcross]');
if isempty(filtfs), filtfs = eqar(1).par_dtstar_specR.filtfs; end

%% Find the average spectrum and use that as a reference
specss_ref = zeros(size(eqar(1).specss));
kk = 0;
for is = 1:length(eqar)
    if isempty(eqar(is).specss), continue; end
    specss_ref = specss_ref + eqar(is).specss;
    kk = kk+1;
end
specss_ref = specss_ref/kk;

%% ===================================================================
%% ------------------ MAP WITH DTSTAR FOR THIS EVENT -----------------
%% ===================================================================

figure(32), clf, hold on
mkfig_CascMAP
set(gca,'xlim',[-134,-120],'ylim',[42,50])
scatter([eqar.slon],[eqar.slat],300,[eqar.dtstar_comb],'filled','MarkerEdgeColor','k')
% % see where the slab depth section is
% scatter(slabdepth(:,3),slabdepth(:,4),10,.05*slabdepth(:,2))

% text([eqar.slon]+0.2,[eqar.slat],{eqar.sta})
cmap = parula;
colormap(cmap)
caxis(tstlim)

% plot section
plot(section_lola(:,1),section_lola(:,2),'k','LineWidth',2)
plot(linterp(section_x,section_lola(:,1),[-100:100:700]'),...
     linterp(section_x,section_lola(:,2),[-100:100:700]'),...
     '.k','MarkerSize',25)
text(linterp(section_x,section_lola(:,1),[-100:100:700]'),...
     linterp(section_x,section_lola(:,2),[-100:100:700]')+0.12,...
     num2str([-100:100:700]'),'FontWeight','bold')

%% colour bar
cbar_custom(gca, 'location',[-133.1 -132.7 43 46],'tickside','right',...
    'lims',tstlim,'tickvals',[tstlim(1):0.5:tstlim(2)],'cmap',cmap,...
    'FontSize',12,'FontWeight','bold',...
	'title',sprintf('$\\Delta t^*_%s$ \\,(s)',eqar(1).phase),'interpreter','latex');

%% ===================================================================
%% -------------------------- DO SOME CALCS --------------------------
%% ===================================================================

%% work out distance to ridge
Xrdg = dist_to_rdg([eqar.slat]',[eqar.slon]'); 
eqar(1).Xrdg = Xrdg(1);
eqar = dealto(eqar,'Xrdg',Xrdg);
[~,Xord] = sort([eqar.Xrdg]); % sort by longitude
[~,Xord] = sort(Xord);
Xrdglims = [min(abs([eqar.Xrdg])),max([eqar.Xrdg])];
kpd = 78; % km per degree of longitude

%% load sed thickness
sectseds = dlmread('mapdata/2D_bathym_section_sediments.txt','\t',2,0); % load sediment thickness for section
sectseds = sectseds(:,[1,2,4]);

%% load + calc. slab depths
slabdepth = dlmread('/Users/zeilon/Work/CASCADIA/plots/depth_to_slab/slabgeom.txt',' ');
% prof_centre=-124.18427/47.09777, az = 100
[slabdepth(:,4),slabdepth(:,3)] = reckon_km(47.09777,-124.18427,slabdepth(:,1),100);
slabdepth(:,1) = slabdepth(:,1)-slabdepth(1,1)+250;
load('mapdata/2D_section_fakeslabdepth.mat'); slbz = faked_slab_depth;

%% Calc differential travel times just from bathymetry + seds
z_seflr = section_z;
sedZ = sectseds(:,3);
sedZ(section_x>320) = 1e3;
z_sedbo = z_seflr - sedZ;
% Make sectseds at least 1km thick on the shelf
% look at flueh
z_slbto = linterp(faked_slab_depth(:,1),faked_slab_depth(:,2),section_x);
z_ocmoh = z_slbto - 6.5e3;
z_fifty = -50e3*ones(size(z_seflr));

% sta_z_seflr = linterp(section_x,z_seflr,Xrdg');
% sta_z_sedbo = linterp(section_x,z_sedbo,Xrdg');
% sta_z_slbto = linterp(section_x,z_slbto,Xrdg');
% sta_z_ocmoh = linterp(section_x,z_ocmoh,Xrdg');
% sta_z_fifty = -50e3*ones(size(z_seflr));

v_asth = 4.4e3;
v_oc = 3.9e3;
v_cr = 3.7e3;
v_sed = 2e3;

TT = (z_seflr-z_sedbo)/v_sed + ...
     (z_sedbo-z_slbto)/v_cr + ...
     (z_slbto-z_ocmoh)/v_oc + ...
     (z_ocmoh-z_fifty)/v_asth;
TT = TT - nanmean(TT) + 1.65;

%% ===================================================================
%% --------------- SECTION WITH DTSTAR FOR THIS EVENT ----------------
%% ===================================================================
figure(17), clf, set(gcf,'pos',[59   258   800   913])
%% topo
maxy = 3000;
miny = -7000;
subplot(3,1,1), hold on
% fill in some layers
fill([section_x;flipud(section_x)],[zeros(size(section_x));flipud(z_seflr)],[221 253 251]/255,'EdgeColor','none')
fill([section_x;flipud(section_x)],[z_seflr;flipud(z_sedbo)],[253 240 221]/255,'EdgeColor','none')
fill([section_x;flipud(section_x)],[z_sedbo;flipud(z_slbto)],[255 241 236]/255,'EdgeColor','none')
fill([section_x;flipud(section_x)],[z_slbto;-1e5*ones(size(section_x))],[245 255 236]/255,'EdgeColor','none')
% draw bathym, seds, slab
plot(slbz(:,1),slbz(:,2),'--k','LineWidth',1)
plot(section_x,z_sedbo,':k','LineWidth',1)
plot(section_x,z_seflr,'k','LineWidth',2)


% figure things
xlim([-80 460]),ylim([miny maxy]+[-1 1])
set(gca,'visible','off','fontsize',12)
% draw axes
plot([-80 720],[0 0],'--k') % x-axis
plot([-80 -80],[miny maxy],'k','LineWidth',2)
plot([-80 -72],[miny miny],'k',[-80 -72],[maxy maxy],'k','LineWidth',1)
text(-80,maxy+1e3,num2str(maxy),'Fontsize',10,'HorizontalAlignment','center','VerticalAlignment','middle')
text(-80,miny-1e3,num2str(miny),'Fontsize',10,'HorizontalAlignment','center','VerticalAlignment','middle')
text(-115,0,'Elev. (m)','Fontsize',17,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','Rotation',90)
% title
text(180,miny-2e3,sprintf('%s $~$  $%s$-wave, %s-comp',epoch2str(evtime,'%Y-%m-%d %H:%M:%S'),phase,comp),'Fontsize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex')
% ridge axis label
plot(0,2500,'vk','MarkerSize',8,'MarkerFaceColor','k')
text(0,3500,'Axis','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14,'Interpreter','latex')
% deformation front label
plot(305,2500,'vk','MarkerSize',8,'MarkerFaceColor','k')
text(305,3500,'DF','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14,'Interpreter','latex')
% Coastline label
plot(440,2500,'vk','MarkerSize',8,'MarkerFaceColor','k')
text(440,3500,'Coastline','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14,'Interpreter','latex')
% Arc label
plot(616,2500,'vk','MarkerSize',8,'MarkerFaceColor','k')
text(616,3500,'Arc','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14,'Interpreter','latex')
% Slab top label
text(330,-6000,'Top of slab','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14,'Interpreter','latex','rotation',-24)
% Stations
plot(Xrdg*kpd,[eqar.selev],'^b','MarkerSize',8,'MarkerFaceColor','r')

% delta tstar
subplot(3,1,2), set(gca,'fontsize',12), hold on
scatter(Xrdg*kpd,[eqar.dtstar],20,'g','MarkerFaceColor','g')
scatter(Xrdg*kpd,[eqar.dtstar_comb],20./[eqar.stds_comb],'b','MarkerFaceColor','b')
% plot(Xrdg(ilan)*kpd,[eqar(ilan).dtstar],'.r','MarkerSize',18)
% figure things
grid on
xlim([-80 460]), ylim([-3 3])
text(-115,0,'$\Delta t^{\ast}_S$','Fontsize',20,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','Rotation',90)

% delta t
subplot(3,1,3), set(gca,'fontsize',12), hold on
plot(section_x,TT,'--r','LineWidth',1.5) % < plot on the prediction of TT from just topog + seds.
scatter(Xrdg*kpd,[eqar.dT],80,'b','MarkerFaceColor','b')
% plot(Xrdg*kpd,[eqar.dT_comb],'.b','MarkerSize',18)
% plot(Xrdg(ilan)*kpd,[eqar(ilan).dT],'.r','MarkerSize',18)
% figure things
grid on
xlim([-80 460]), ylim([-2.5 2.5])
ax = axis;
text(-115,0,'$\delta T_S$','Fontsize',20,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','Rotation',90)
text(mean(ax(1:2)),-3.5,'Distance from ridge (km)','Fontsize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex')

%% ===================================================================
%% ---------------------------- WAVEFORMS ----------------------------
%% ===================================================================

for is = 1:length(eqar)
    dat = eqar(is).(['dat',comp])';

    W = 2*filtfs./eqar(is).samprate;
    [fb, fa] = butter(2,W);
    dd(:,is) = filtfilt(fb,fa,dat);
end

%    plot
figure(2), clf, set(gcf,'position',[150 000 600 800]), hold on
for is = 1:length(eqar)
	tt = eqar(is).tt-eqar(is).abs_arrT;

    hp = plot(tt,amp_seis*dd(:,is)/max(max(abs(dd))) + Xord(is),'LineWidth',2);
    set(hp,'color',colour_get(abs(eqar(is).Xrdg),Xrdglims(2),Xrdglims(1),flipud(jet)))
    text(datwind(1)+2,Xord(is)+0.4,eqar(is).sta,...
        'FontSize',8,'interpreter','latex','HorizontalAlignment','left')
    text(datwind(2) + 1,Xord(is),sprintf('%.0f km',eqar(is).Xrdg*kpd),...
        'FontSize',10,'interpreter','latex')
end
plot([1;1]*window,[0;length(eqar)+1]*[1 1],'--k','LineWidth',1)
set(gca,'FontSize',12,'YTick',[],'position',[0.16 0.11 0.75 0.815])
axis([datwind,-0.5,length(eqar)+1.5])
xlabel('Time from phase onset (s)','FontSize',14,'interpreter','latex')
title(sprintf('%s $~$ %s-wave, %s-comp, $f_{hi}~$: %.2f',epoch2str(evtime,'%Y-%m-%d %H:%M:%S'),...
    phase,comp,fmax),'Fontsize',16,'Interpreter','latex')

% ridge axis label
text(datwind(1)-2,sum([eqar.Xrdg] < 0)+0.5,'\textbf{Axis $\succ$}',...
    'HorizontalAlignment','right','VerticalAlignment','bottom',...
    'FontSize',15,'FontWeight','bold','Interpreter','latex')
% deformation front label
text(datwind(1)-2,sum([eqar.Xrdg] < 305/kpd)+0.5,'\textbf{DF $\succ$}',...
    'HorizontalAlignment','right','VerticalAlignment','bottom',...
    'FontSize',15,'FontWeight','bold','Interpreter','latex')
% Coastline label
text(datwind(1)-2,sum([eqar.Xrdg] < 440/kpd)+0.5,'\textbf{Coast $\succ$}',...
    'HorizontalAlignment','right','VerticalAlignment','bottom',...
    'FontSize',15,'FontWeight','bold','Interpreter','latex')
% Arc label
% text(-32,sum([eqar_obs.Xrdg] < 616/kpd)+0.5,'Arc $\rightarrow$',...
%     'HorizontalAlignment','right','VerticalAlignment','bottom',...
%     'FontSize',15,'FontWeight','bold','Interpreter','latex')


% %% ===================================================================
% %% ------------------------ SPECTRA (W/ FITS) ------------------------
% %% ===================================================================
% figure(3), clf, set(gcf,'position',[200 10 800 600]), hold on
% for is = 1:length(eqar_plot)
%     
%     frq = eqar_plot(is).frq;
%     lnR = log(eqar_plot(is).specss./specss_ref);
% 
%     ind = frq <= par.hifrq & frq >= par.lofrq;
% 
%     fo = fit(frq(ind),lnR(ind),'poly1');
%   
%     hp = plot(frq(ind),lnR(ind) - fo.p2,'o-','LineWidth',1.5);
%     hpf = plot(frq(ind),fo.p1*frq(ind),'--');
% 
%     set([hp,hpf],'color',colour_get(abs(eqar_plot(is).Xrdg),Xrdglims_plot(2),Xrdglims_plot(1),flipud(jet)))
%     xlim([0 par.hifrq])
%     
% end
% xlim([par.lofrq par.hifrq+0.01])
% % ylim([-2.5 3])
% set(gca,'FontSize',8)
% xlabel('Frequency (Hz)','interpreter','latex','FontSize',16)
% ylabel('Spectral ratio: $~\ln~(A_i/A_0)$','interpreter','latex','FontSize',16)
% title('OBS spectral ratios relative to mean spectrum','interpreter','latex','FontSize',18)
% legend off
% 
% cbar_custom(gca, 'location',[0.06 .12 1 1.18],'tickside','bottom',...
%     'lims',kpd*Xrdglims_plot,'tickvals',round_level([kpd*Xrdglims_plot(1):50:kpd*Xrdglims_plot(2)],100),...
%     'FontSize',10,'FontWeight','bold','cmap',flipud(jet),...
%     'title','Distance E of ridge (km)','interpreter','latex');
% 



end

