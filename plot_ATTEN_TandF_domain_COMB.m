function plot_ATTEN_TandF_domain_COMB( eqar, plotopt )
% plotopt can be "all" or "obs"

if nargin<2
    plotopt = 'all';
end

datwind = [-70 100];
filtfs = [0.05 2]; % in hz
tstlim = 1.5*[-1 1];
% refsta = 'LON';

t = tauptime('p',eqar(1).phase,'deg',eqar(1).gcarc); 
evtime = eqar(1).pred_arrT - t(1).time;

% get "good" stations
gd=[]; for is = 1:length(eqar), if ~isempty(eqar(is).dtstar_comb) && ~isnan(eqar(is).dtstar_comb), gd = [gd;is]; end, end
eqar_gd = eqar(gd);
% put NaNs in for stations without specR dtstar:
for is=1:length(eqar_gd), if isempty(eqar_gd(is).dtstar), eqar_gd(is).dtstar=NaN; end; end

%% get deets
phase = eqar_gd(1).phase;
comp = eqar_gd(1).par_dtstar.comp;
window = [eqar_gd(1).par_dtstar_comb.wind.prex,eqar_gd(1).par_dtstar_comb.wind.postx];
fmax = eqar_gd(1).par_dtstar_comb.inv.fmax;
if isempty(filtfs), filtfs = eqar_gd(1).par_dtstar_specR.filtfs; end

%% parse obs/land stations
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
fprintf('Using %s database\n',[dbdir,dbnam])
db = dbopen([dbdir,dbnam],'r');
dbsi = dblookup_table(db,'site');
[stas,statypes] = dbgetv(dbsi,'sta','statype');
dbclose(db);
[~,inds,~] = intersect(stas,{eqar_gd.sta},'stable');
ilan = strcmp(statypes(inds),'LAND');
iobs = strcmp(statypes(inds),'OBS');

% % << MAKE DEFUNCT>>
% iobs = ~cellfun('isempty',regexp({eqar_gd.sta},'([J,F,M,G]*[0-9][0-9][A-C]$)'));
% ilan =  cellfun('isempty',regexp({eqar_gd.sta},'([J,F,M,G]*[0-9][0-9][A-C]$)'));
% 
% % account for stupidly named land stations
% ilan(find(strcmp({eqar_gd(iobs).sta},'M02C') & [eqar_gd(iobs).slat]==41.392))=true;
% iobs(find(strcmp({eqar_gd(iobs).sta},'M02C') & [eqar_gd(iobs).slat]==41.392))=false;


% if nargin<2 || isempty(refsta)
% mlat = mean([eqar_gd.slat]); 
% mlon = mean([eqar_gd.slon]); 
% dlalo = distance(mlat,mlon,[eqar_gd.slat],[eqar_gd.slon]);
% iref = find(dlalo==min(dlalo));
% else
% iref = find(strcmp({eqar_gd.sta},refsta));
% end

%% Find the average spectrum and use that as a reference
specss_ref = zeros(size(eqar_gd(1).specss));
kk = 0;
for is = 1:length(eqar_gd)
    if isempty(eqar_gd(is).specss), continue; end
    specss_ref = specss_ref + eqar_gd(is).specss;
    kk = kk+1;
end
specss_ref = specss_ref/kk;


%% work out distance to ridge
roughjdf = [-130.5,44.36;-128.7,48.98];
Xrdg = dist2line(roughjdf(1,:),roughjdf(2,:),[[eqar_gd.slon]',[eqar_gd.slat]'])'; 
eqar_gd(1).Xrdg = Xrdg(1);
eqar_gd = dealto(eqar_gd,'Xrdg',Xrdg);
[~,Xord] = sort([eqar_gd.Xrdg]); % sort by longitude
[~,Xord] = sort(Xord);
Xrdglims = [min(abs([eqar_gd.Xrdg])),max([eqar_gd.Xrdg])];
kpd = 78; % km per degree of longitude

eqar_obs = eqar_gd(iobs);
[~,Xord_obs] = sort([eqar_obs.Xrdg]); % sort by longitude
[~,Xord_obs] = sort(Xord_obs);
Xrdglims_obs = [min(abs([eqar_obs.Xrdg])),max([eqar_obs.Xrdg])];
Xrdglims_obs = [0,max([eqar_obs.Xrdg])];

if strcmp(plotopt,'obs')
    Xord_plot = Xord_obs;
    eqar_plot = eqar_obs;
    Xrdglims_plot = Xrdglims_obs;
elseif strcmp(plotopt,'all')
    Xord_plot = Xord;
    eqar_plot = eqar_gd;
    Xrdglims_plot = Xrdglims;    
end
 
%% ===================================================================
%% ------------------ MAP WITH DTSTAR FOR THIS EVENT -----------------
%% ===================================================================

figure(32), clf, hold on
mkfig_CascMAP
scatter([eqar_gd.slon],[eqar_gd.slat],200,[eqar_gd.dtstar_comb],'filled')
% text([eqar_gd.slon]+0.2,[eqar_gd.slat],{eqar_gd.sta})
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
cbar_custom(gca, 'location',[-132.4 -132 43 46],'tickside','right',...
    'lims',tstlim,'tickvals',[tstlim(1):0.5:tstlim(2)],'cmap',cmap,...
    'FontSize',12,'FontWeight','bold',...
	'title',sprintf('$\\Delta t^*_%s$ \\,(s)',eqar(1).phase),'interpreter','latex');
 

%% ===================================================================
%% --------------- SECTION WITH DTSTAR FOR THIS EVENT ----------------
%% ===================================================================
figure(17), clf, set(gcf,'pos',[59   258   871   613])
%% topo
subplot(5,1,1), hold on
plot(section_x,section_z,'k','LineWidth',2)
% figure things
xlim([-80 720]),ylim([-5001 5001])
set(gca,'visible','off','fontsize',12)
% draw axes
plot([-80 720],[0 0],'--k') % x-axis
plot([-80 -80],[-5000 5000],'k','LineWidth',2)
plot([-80 -72],[-5000 -5000],'k',[-80 -72],[5000 5000],'k','LineWidth',1)
text(-80,7000,'5000','Fontsize',10,'HorizontalAlignment','center','VerticalAlignment','middle')
text(-80,-7000,'-5000','Fontsize',10,'HorizontalAlignment','center','VerticalAlignment','middle')
text(-115,0,'Elev (m)','Fontsize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','Rotation',90)
% title
text(300,-7000,sprintf('%s $~$  $%s$-wave, %s-comp',epoch2str(evtime,'%Y-%m-%d %H:%M:%S'),phase,comp),'Fontsize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex')
% ridge axis label
plot(0,4500,'vk','MarkerSize',8,'MarkerFaceColor','k')
text(0,6500,'Axis','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'Interpreter','latex')
% deformation front label
plot(305,4500,'vk','MarkerSize',8,'MarkerFaceColor','k')
text(305,6500,'DF','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'Interpreter','latex')
% Coastline label
plot(440,4500,'vk','MarkerSize',8,'MarkerFaceColor','k')
text(440,6500,'Coastline','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'Interpreter','latex')
% Arc label
plot(616,4500,'vk','MarkerSize',8,'MarkerFaceColor','k')
text(616,6500,'Arc','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'Interpreter','latex')


% delta tstar
subplot(5,1,2:3), set(gca,'fontsize',12), hold on
plot(Xrdg(iobs)*kpd,[eqar_gd(iobs).dtstar],'.g','MarkerSize',10)
plot(Xrdg(iobs)*kpd,[eqar_gd(iobs).dtstar_comb],'.b','MarkerSize',18)
plot(([eqar_gd(ilan).slon]+124)*kpd + 440,[eqar_gd(ilan).dtstar],'.g','MarkerSize',10)
plot(([eqar_gd(ilan).slon]+124)*kpd + 440,[eqar_gd(ilan).dtstar_comb],'.r','MarkerSize',18)
% plot(Xrdg(ilan)*kpd,[eqar_gd(ilan).dtstar],'.r','MarkerSize',18)
% figure things
grid on
text(-115,0.5,'$\Delta t^{\ast}_S$','Fontsize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','Rotation',90)
xlim([-80 720])

% delta t
subplot(5,1,4:5), set(gca,'fontsize',12), hold on
plot(Xrdg(iobs)*kpd,[eqar_gd(iobs).dT],'.g','MarkerSize',10)
plot(Xrdg(iobs)*kpd,[eqar_gd(iobs).dT_comb],'.b','MarkerSize',18)
plot(([eqar_gd(ilan).slon]+124)*kpd + 440,[eqar_gd(ilan).dT],'.g','MarkerSize',10)
plot(([eqar_gd(ilan).slon]+124)*kpd + 440,[eqar_gd(ilan).dT_comb],'.r','MarkerSize',18)
% plot(Xrdg(ilan)*kpd,[eqar_gd(ilan).dT],'.r','MarkerSize',18)
% figure things
grid on
text(-115,0,'$\Delta t_S$','Fontsize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','Rotation',90)
text(320,-4.5,'Distance from ridge (km)','Fontsize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex')
xlim([-80 720]), ylim([-3 3])

%% ===================================================================
%% ---------------------------- WAVEFORMS ----------------------------
%% ===================================================================

for is = 1:length(eqar_plot)
    dat = eqar_plot(is).(['dat',comp])';

    W = 2*filtfs./eqar_plot(is).samprate;
    [fb, fa] = butter(2,W);
    dd(:,is) = filtfilt(fb,fa,dat);
end

%    plot
figure(2), clf, set(gcf,'position',[150 000 600 800]), hold on
for is = 1:length(eqar_plot)
	tt = eqar_plot(is).tt-eqar_plot(is).abs_arrT;

    hp = plot(tt,5*dd(:,is)/max(max(abs(dd))) + Xord_plot(is),'LineWidth',2);
    set(hp,'color',colour_get(abs(eqar_plot(is).Xrdg),Xrdglims_plot(2),Xrdglims_plot(1),flipud(jet)))
    text(datwind(1)+2,Xord_plot(is)+0.4,eqar_plot(is).sta,...
        'FontSize',8,'interpreter','latex','HorizontalAlignment','left')
    text(datwind(2) + 1,Xord_plot(is),sprintf('%.0f km',eqar_plot(is).Xrdg*kpd),...
        'FontSize',10,'interpreter','latex')
end
plot([1;1]*window,[0;length(eqar_plot)+1]*[1 1],'--k','LineWidth',1)
set(gca,'FontSize',12,'YTick',[],'position',[0.16 0.11 0.75 0.815])
axis([datwind,-0.5,length(eqar_plot)+1.5])
xlabel('Time from phase onset (s)','FontSize',14,'interpreter','latex')
title(sprintf('%s $~$ %s-wave, %s-comp, $f_{hi}~$: %.2f',epoch2str(evtime,'%Y-%m-%d %H:%M:%S'),...
    phase,comp,fmax),'Fontsize',16,'Interpreter','latex')

% ridge axis label
text(datwind(1)-2,sum([eqar_plot.Xrdg] < 0)+0.5,'\textbf{Axis $\succ$}',...
    'HorizontalAlignment','right','VerticalAlignment','bottom',...
    'FontSize',15,'FontWeight','bold','Interpreter','latex')
% deformation front label
text(datwind(1)-2,sum([eqar_plot.Xrdg] < 305/kpd)+0.5,'\textbf{DF $\succ$}',...
    'HorizontalAlignment','right','VerticalAlignment','bottom',...
    'FontSize',15,'FontWeight','bold','Interpreter','latex')
% Coastline label
text(datwind(1)-2,sum([eqar_plot.Xrdg] < 440/kpd)+0.5,'\textbf{Coast $\succ$}',...
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

