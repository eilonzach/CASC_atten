function plot_ATTEN_TandF_domain( eqar,evtime )

datwind = [-30 50];
filtfs = [0.1 10];
refsta = 'LON';

% get "good" stations
gd=[]; for is = 1:length(eqar), if ~isempty(eqar(is).dtstar), gd = [gd;is]; end, end

iref = find(strcmp({eqar.sta},refsta));

%% work out distance to ridge
roughjdf = [-130.5,44.36;-128.7,48.98];
Xrdg = dist2line(roughjdf(1,:),roughjdf(2,:),[[eqar.slon]',[eqar.slat]'])'; 
eqar(1).Xrdg = Xrdg(1);
eqar = dealto(eqar,'Xrdg',Xrdg);
[~,Xord] = sort([eqar(gd).Xrdg]); % sort by longitude
[~,Xord] = sort(Xord);

%% ------------------ MAP WITH DTSTAR FOR THIS EVENT ------------------

figure(31), clf, hold on
mkfig_CascMAP
scatter([eqar(gd).slon],[eqar(gd).slat],180,[eqar(gd).dtstar],'filled')
% scatter([eqar(gd).slon],[eqar(gd).slat],90,-[eqar(gd).Xrdg],'filled')
text([eqar(gd).slon]+0.2,[eqar(gd).slat],{eqar(gd).sta})
colorbar


par = eqar(gd(1)).par_dtstar;
if isempty(filtfs), filtfs = par.filtfs; end



%% ------------------ WAVEFORMS AND SPECTRA (W/ FITS) ------------------
figure(2)
clf
set(gcf,'position',[150 000 800 800])

%% WAVEFORMS
subplot(211), hold on
for ig = 1:length(gd)
    is = gd(ig);

    tt = eqar(is).tt-eqar(is).abs_arrT;
    dat = eqar(is).(['dat',par.comp])';

    W = 2*filtfs./eqar(is).samprate;
    [fb, fa] = butter(2,W);
    dd = filtfilt(fb,fa,dat);

    hp = plot(tt,dd/max(abs(dd)) + Xord(ig),'LineWidth',2);
    set(hp,'color',colour_get(eqar(is).Xrdg,min([eqar(gd).Xrdg]),max([eqar(gd).Xrdg])))
    text(datwind(1)-2,Xord(ig),eqar(is).sta,'FontSize',11,'interpreter','latex','HorizontalAlignment','right')
    text(datwind(2) + 1,Xord(ig),[num2str(eqar(is).Xrdg,2) '$^{\circ}$'],'FontSize',10,'interpreter','latex')
end
plot([1;1]*par.window,[0;length(gd)+1]*[1 1],'--k','LineWidth',1)
set(gca,'FontSize',8,'YTick',[])
axis([datwind,-0.5,length(gd)+1.5])
xlabel('Time from phase onset (s)','FontSize',14,'interpreter','latex')
title(sprintf('%s $~$ Comp: %s, $f_{hi}~$: %.1f',epoch2str(evtime,'%Y-%m-%d %H:%M:%S'),par.comp,par.hifrq),'Fontsize',14,'Interpreter','latex')


specss_ref = eqar(iref).specss;
           
%% SPECTRA + FITS
subplot(212), hold on
for ig = 1:length(gd)
    is = gd(ig);
    
    frq = eqar(is).frq;
    lnR = log(eqar(is).specss./specss_ref);
    
    ind = frq < par.hifrq;
    hp = plot(frq(ind),lnR(ind),'o-','LineWidth',1.5);
    hpf = plot(fit(frq(ind),lnR(ind),'poly1'),'--');

    set([hp,hpf],'color',colour_get(eqar(is).Xrdg,min([eqar(gd).Xrdg]),max([eqar(gd).Xrdg])))
    xlim([0 par.hifrq])
    
end
xlim([0.03 par.hifrq+0.01])
ylim([-2 3])
set(gca,'FontSize',8)
xlabel('Frequency (Hz)','interpreter','latex','FontSize',14)
ylabel('Spectral ratio: $~\ln~(A_i/A_0)$','interpreter','latex','FontSize',14)
legend off



end

