clear all

refsta = 'RADR';
component = 'T'; %'Z', 'R', or 'T'
snrmin = 0;
wlen = 15;
t0 = -5;
hifrq = 0.2; % uppermost freq to fit (Hz)
mavwind = 1; % length of moving average to smooth spectrum (1 for no smoothing)
latlim = [42.5 48];
lonlim = [-131 -120];

load eg_cleandata.mat
load('coast')
jdf = dlmread('mapdata/ridge_xy');
fzs = dlmread('mapdata/transforms_xy');

%% only use data from swath of stations
gd = [data.slat] < latlim(2) & [data.slat] > latlim(1) & [data.slon] < lonlim(2) & [data.slon] > lonlim(1);
data = data(gd);
eval(sprintf('for is = 1:length(data), data(is).dat = data(is).data%s; end',component))

figure(9), clf, hold on
set(gcf,'pos',[ 10  50  900  780])
scatter([data.slon],[data.slat],70,colour_get([data.slon],max([data.slon]),min([data.slon])),'filled')
text([data.slon],[data.slat]+0.1,{data.stnm})
plot(long,lat)
scatter( -121.760374, 46.852886,150,'k','^','filled') % Mt Ranier
plot(jdf(:,1),jdf(:,2),'k','Linewidth',2)
plot(fzs(:,1),fzs(:,2),'--k','Linewidth',1.5)
xlim(lonlim + [-0.5 0.5])
ylim(latlim + [-1 1])

% Notes from geoff's quickq:
%  Need at least 5 Hz samprate
%  Sample over 4s window
%  begin data window so phase arrival is 10% into window, and noise window
%   is the 4s before that... but NB this wasn't for teleseisms



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

%% take spectral ratio relative to reference station
for is = 1:nstas
    data(is).lnR = log(data(is).specss./data(irst).specss); % << here using smoothed spectrum 
end
    
figure(2), clf, hold on
for is = 1:nstas
if ~any(is==gd), continue, end    
ind = frq < hifrq;
hp = plot(frq(ind),data(is).lnR(ind),'o-');
set(hp,'color',colour_get(data(is).slon,max([data.slon]),min([data.slon])))
xlim([0 hifrq])
fo = fit(frq(ind),data(is).lnR(ind),'poly1');
data(is).dtstar = fo.p1/(-pi);
data(is).dtstar_cb = confint(fo); data(is).dtstar_cb = data(is).dtstar_cb(1,:)'; % shouldn't this have a /(pi);?
end

figure(3), clf, hold on
plot([data(gd).slon],[data(gd).dtstar],'.','MarkerSize',30)
text([data(gd).slon],[data(gd).dtstar]+0.02,{data(gd).stnm})
eb = [data(gd).dtstar_cb]';
errorbar([data(gd).slon]',[data(gd).dtstar]',eb(:,1),eb(:,2),'o')
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
set(hps,'color',colour_get(data(is).slon,max([data.slon]),min([data.slon])))
% xlim(ax(1),[0 0.8*nyq]);
xlabel('Hz'); 
ylabel('amplitude, nm/Hz')
% a=axis(ax(1));
% % ylim_max=10^(floor(log10(max(specs(1:n4))))+1);
% % ylim(ax(1),ylim_max.*[1.e-4 1])
xlim([0 hifrq]);
ylim([-10,-3]);

end

