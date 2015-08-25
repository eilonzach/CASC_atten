% plot waveforms to see atten. in waveform domain
filtfreqs = [40 0.1];

%% select stations we want
badstas = {'M03C','M05C','J25C'};
gdgd = intersect(gd,find([data.elev]<0));
for ii = 1:length(badstas)
    gdgd(strcmp({data(gdgd).stnm},badstas{ii}))=[];
end

%% Work out some event data
[gcarc,baz] = distance(mean([data.slat]),mean([data.slon]),evinfo.Origins(1).Latitude,evinfo.Origins(1).Longitude);




figure(19),clf,
set(gcf,'position',[0 0 800 800])
jj = find(tt >= t0-2*wlen & tt < t0 + 3*wlen);
re_dt = unique(round_level(diff(tt),0.0001));
W = 2*re_dt./filtfreqs;
[fb, fa] = butter(2,W);

[~,Xord] = sort([data(gdgd).Xrdg]); % sort by longitude
[~,Xord] = sort(Xord);

kk = 0;
subplot(211), hold on
for is = 1:nstas
if ~any(is==gdgd), continue, end    
kk = kk+1;
dd = filtfilt(fb,fa,data(is).dat);

hp = plot(tt(jj),dd(jj)/max(abs(dd(jj))) + Xord(kk),'LineWidth',2);
set(hp,'color',colour_get(data(is).Xrdg,min([data(gdgd).Xrdg]),max([data(gdgd).Xrdg])))
text(-wlen-2,Xord(kk),data(is).stnm,'FontSize',11,'interpreter','latex','HorizontalAlignment','right')
text(2*wlen + 1,Xord(kk),[num2str(data(is).Xrdg,2) '$^{\circ}$'],'FontSize',10,'interpreter','latex')
xlim([-wlen wlen*2])
end
plot([1;1]*[t0 t0+wlen],[0;kk+1]*[1 1],'--k','LineWidth',1)
set(gca,'FontSize',8,'YTick',[])
xlabel('Time from phase onset (s)','FontSize',11,'interpreter','latex')
title(sprintf('%s $~$ Comp: %s, $f_{hi}~$: %.1f',evinfo.PreferredTime(1:19),component,hifrq),'Fontsize',14,'Interpreter','latex')

            
subplot(212), hold on
for is = 1:nstas
if ~any(is==gdgd), continue, end    
ind = frq < hifrq;
hp = plot(frq(ind),data(is).lnR(ind),'o-','LineWidth',1.5);
set(hp,'color',colour_get(data(is).Xrdg,min([data(gdgd).Xrdg]),max([data(gdgd).Xrdg])))
xlim([0 hifrq])
hp = plot(fit(frq(ind),data(is).lnR(ind),'poly1'),'--');
set(hp,'color',colour_get(data(is).Xrdg,min([data(gdgd).Xrdg]),max([data(gdgd).Xrdg])))
end
xlim([0.03 hifrq+0.01])
ylim([-1 1.5])
set(gca,'FontSize',8)
xlabel('Frequency (Hz)','interpreter','latex','FontSize',11)
ylabel('Spectral ratio: $~\ln~(A_i/A_0)$','interpreter','latex','FontSize',11)
% xlabel
legend off
            
save2pdf(19,'OBS_waves_spectra','figs');