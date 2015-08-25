clear all

LAND = {'FORK','WISH','D03D','DOSE','GNW','D04E','SP2','STOR','RATT','OBSR','LON','LTY'};
   
refsta = LAND{1};

component = 'T'; %'Z', 'R', or 'T'
snrmin = 0;
wlen = 40;
t0 = -20;
hifrq = 0.1; % uppermost freq to fit (Hz)
mavwind = 1; % length of moving average to smooth spectrum (1 for no smoothing)
latlim = [42 50];
lonlim = [-131 -120];

load eg_cleandata.mat
jdf = dlmread('mapdata/ridge_xy');

%% Work out some event data
load eg_cleaninfo.mat
evtime = evinfo.Origins(1).Time(1:10);
[gcarc,baz] = distance(mean([data.slat]),mean([data.slon]),evinfo.Origins(1).Latitude,evinfo.Origins(1).Longitude);


%% only use data from OBS
[~,gd] = intersect({data.stnm},LAND);
data = data(gd);

%% work out distance to ridge
for id = 1:length(data)
    data(id).Xrdg = min(distance(data(id).slat,data(id).slon,jdf(:,2),jdf(:,1)));
end


figure(9), clf, hold on
set(gcf,'pos',[ 10  50  900  780])
CascMAP
scatter([data.slon],[data.slat],90,colour_get([data.Xrdg],min([data.Xrdg]),max([data.Xrdg])),'filled')
text([data.slon]+0.1,[data.slat]+0.1,{data.stnm},'FontWeight','bold')

samprate = 1./unique(round_level(diff(data(1).timeaxis),0.001));
nstas = length(data);
dt = 1./samprate;
nyq=0.5./dt;

tt = data(1).timeaxis;
jn = find(tt >= t0-wlen & tt < t0);
js = find(tt >= t0 & tt < t0 + wlen);
nw=length(jn);
wdo=tukeywin(nw,0.2)';   % 0.2 gives 10% cosine taper

%% Reference station
irst = find(strcmp(refsta,{data.stnm}));
plot(data(irst).slon, data(irst).slat,'ok','MarkerSize',10,'LineWidth',2.5)

figure(1), clf
for is = 1:nstas
dat = data(is).datfilt;
dn = dat(jn);
ds = dat(js);
data(is).snr = var(ds)/var(dn);    % signal-to-noise, useful for filtering fits
% detrend
dn = detrend(dn);
ds = detrend(ds);
% window
dn=dn.*wdo;
ds=ds.*wdo;
% pad
% dn = [zeros(1,10),dn,zeros(1,10)];
% ds = [zeros(1,10),ds,zeros(1,10)];

ntap=2;
nft=2^nextpow2(nw);
[specn,~]=pmtm(dn,ntap,nft,1./dt);
[specs,frq]=pmtm(ds,ntap,nft,1./dt);
frq=frq(2:length(frq));
specn=specn(2:length(specn));
specs=specs(2:length(specs));

%convert power to amplitude, and integrate to displacement
data(is).specn=(specn.^0.5)./(2.*pi.*frq);
data(is).specs=(specs.^0.5)./(2.*pi.*frq);

data(is).specss = moving_average(data(is).specs,mavwind);

data(is).fmax = frq(max([find(data(is).specs<data(is).specn,1,'first'),2])-1);
if isempty(data(is).fmax), data(is).fmax = 0; end

if data(is).snr < snrmin, continue, end
%% plot data series
figure(1)
subplot(212), hold on
ja=max(js) + (1:nw);
da = detrend(data(is).dat(ja)).*wdo;
hold on;
plot(tt(jn),dn,'r',tt(js),ds,'g',tt(ja),da,'k');
xlabel('time from pick, s')
xlim([min(tt(jn)) max(tt(ja))])

end

%% look at maximum frequencies above noise
figure(87), clf, hold on
fmaxs = [data.fmax]';
hist(fmaxs,100); xlim([0 2])
line(hifrq*[1 1],[0 max(hist(fmaxs))],'Color','r','LineStyle','--','LineWidth',2)
title('maximum freq where signal is above noise')

%% only use good stations stations
gd = find([data.snr]>=snrmin & [data.fmax] >= hifrq);
fprintf('Only %.0f stations left once snrmin and fmax accounted for\n',length(gd))

%% take spectral ratio relative to reference station
for is = 1:nstas
    data(is).lnR = log(data(is).specss./data(irst).specss); % << here using smoothed spectrum 
end
    
figure(2), clf, hold on
for is = 1:nstas
if ~any(is==gd), continue, end    
ind = frq < hifrq;
% plot spectral ratios
hp = plot(frq(ind),data(is).lnR(ind),'o-');
set(hp,'color',colour_get(data(is).Xrdg,min([data.Xrdg]),max([data.Xrdg])),'LineWidth',2)
xlim([0 hifrq])

fo = fit(frq(ind),data(is).lnR(ind),'poly1');
data(is).dtstar = fo.p1/(-pi);
data(is).dtstar_cb = confint(fo); data(is).dtstar_cb = data(is).dtstar_cb(1,:)';
hp = plot(fo);
set(hp,'color',colour_get(data(is).Xrdg,min([data.Xrdg]),max([data.Xrdg])),'LineStyle',':','LineWidth',1.5)

set(legend,'visible','off')

title('Spectral ratios to reference station')
ylabel('Logarithm of spectral ratio')
xlabel('Frequency')
end


figure(3), clf, hold on
plot([data(gd).Xrdg],[data(gd).dtstar],'.','MarkerSize',30)
text([data(gd).Xrdg],[data(gd).dtstar]+0.02,{data(gd).stnm})
eb = [data(gd).dtstar_cb]';
errorbar([data(gd).Xrdg]',[data(gd).dtstar]',eb(:,1),eb(:,2),'o')
grid on
ylabel('$\Delta t^{\ast}$ where negative means more high-f energy','interpreter','latex','Fontsize',15)
    
figure(1)
for is = 1:nstas
if ~any(is==gd), continue, end    
%% spectral plot
subplot(211), hold on
n4=length(find(frq<=0.8*nyq));
hps = plot(frq(1:n4),log10(data(is).specss(1:n4)),'b','LineWidth',1.5); % << here using smoothed spectrum 
hpn = plot(frq(1:n4),log10(data(is).specn(1:n4)),'k','LineWidth',1.5);
set(hps,'color',colour_get(data(is).Xrdg,min([data.Xrdg]),max([data.Xrdg])))
% xlim(ax(1),[0 0.8*nyq]);
xlabel('Hz'); 
ylabel('amplitude, nm/Hz')
% a=axis(ax(1));
% % ylim_max=10^(floor(log10(max(specs(1:n4))))+1);
% % ylim(ax(1),ylim_max.*[1.e-4 1])
xlim([0 hifrq]);
ylim([-10,-3]);

end



% plot waveforms to see atten. in waveform domain
filtfreqs = [40 0.1];

figure(91),clf, hold on
jj = find(tt >= t0-2*wlen & tt < t0 + 3*wlen);
re_dt = unique(round_level(diff(tt),0.0001));
W = 2*re_dt./filtfreqs;
[fb, fa] = butter(2,W);

[~,Xord] = sort([data.Xrdg]); % sort by longitude
[~,Xord] = sort(Xord);

for is = 1:nstas
if ~any(is==gd), continue, end    
dd = filtfilt(fb,fa,data(is).dat);

hp = plot(tt(jj),dd(jj)/max(abs(dd(jj))) + Xord(is),'LineWidth',2);
set(hp,'color',colour_get(data(is).Xrdg,min([data.Xrdg]),max([data.Xrdg])))
text(-wlen-5,Xord(is),data(is).stnm,'FontSize',16)
text(2*wlen + 1,Xord(is),[num2str(data(is).Xrdg,2) '  deg'],'FontSize',14)
xlim([-wlen wlen*2])
end
plot([1;1]*[t0 t0+wlen],[0;length(gd)+1]*[1 1],'--k','LineWidth',1)
set(gca,'FontSize',14,'YTick',[])
xlabel('Time (s) from phase onset','FontSize',15)
title(sprintf('Event %s gcarc=%.1f baz=%.0f filter %.2f to %.2f sec',...
                evtime,gcarc,baz,filtfreqs(1),filtfreqs(2)),'FontSize',18)
