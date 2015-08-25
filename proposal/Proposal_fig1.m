
redo = 1;
if redo
clear all
close all

event = '20140513'; %WAS <==
% event = '20131025'; %< gives massive signal (use 25 s wlen)
phase = 'S';

refsta = 'RADR'; % was RADR
component = 'T'; %'Z', 'R', or 'T'

% parms for data treatment
resamprate = 40 ; % new, common sample rate
filtfs = [25 0.5];
xcorwind = [-10 30];
xcorlagmax = 2.5;
acormin = 0.7;
snrmin = 3;

rmresp = 1;

% parms for spectra calc
% good answer with: wlen = 18; t0 = 0; hifrq = 0.2;
% also reasonable with: wlen = 25; t0 = 0; hifrq = 0.3;
wlen = 18;
t0 = 0;
hifrq = 0.2; % uppermost freq to fit (Hz)
mavwind = 1; % length of moving average to smooth spectrum (1 for no smoothing)
latlim = [43.5 49];
lonlim = [-131 -120];



%% ================================================================
%% ============== shouldn't need to alter below here ==============

addpath /Users/zeilon/Documents/MATLAB/jingle/recreadplot

load(sprintf('eg_evt_sta_info_%s_%s.mat',event,phase));
load(sprintf('eg_traces_%s_%s.mat',event,phase));

evtime = datenum(evinfo.PreferredTime,'yyyy-mm-dd HH:MM:SS');
re_dt = 1./resamprate;

%% remove instrument response, rotate channels, quality control
kk = 0;
for is = 1:size(traces,1)
    if isempty(traces(is,1).network), continue, end
for ic = 1:3
	trace = traces(is,ic);
	% remove instrument response
    if rmresp
    	trace = rm_resp(trace); 
    else
        trace.data_cor = trace.data./max(trace.data);
    end
	
    % resample waveform to ensure same data length
	b = (trace.startTime - evtime )*24*3600  - stinfo(is).Ptime;
	e = (trace.endTime - evtime )*24*3600  - stinfo(is).Ptime;
    dt = (e-b)/trace.sampleCount;
    old_timeaxis = [b:dt:(b + dt*(trace.sampleCount-1))]';
    try
	new_timeaxis = [-stinfo(is).twin(1):re_dt:stinfo(is).twin(2)]';
    catch
    new_timeaxis = [-150:re_dt:150]';
    end

    dat(:,ic) = interp1(old_timeaxis,trace.data_cor,new_timeaxis);	
end

% kill nans
if any(any(isnan(dat))), disp('NANS!'), continue, end

%  rotate to RTZ
    az1 = traces(is,1).azimuth;
    az2 = traces(is,2).azimuth;
        if isempty(az1)
            if     strcmp(traces(is,1).channel(3),'E') az1=90; az2=0;
            elseif strcmp(traces(is,1).channel(3),'N') az1=0;  az2=90;
            else  continue
            end
        end
    baz = azimuth(trace.latitude,trace.longitude,evinfo.PreferredLatitude,evinfo.PreferredLongitude);
	R_az = baz+180;
	T_az = R_az+90;
	dataR = dat(:,1)*cosd(R_az-az1(1))+dat(:,2)*cosd(R_az-az2(1)); % datN*cos(Raz) + datE*cos(Raz)
	dataT = dat(:,1)*cosd(T_az-az1(1))+dat(:,2)*cosd(T_az-az2(1)); % datN*cos(Taz) + datE*cos(Taz)
    dataZ = dat(:,3);
    % filter
    W = 2*re_dt./filtfs;
    [fb, fa] = butter(2,W);

% build up structure 
    kk = kk + 1;
	data(kk).slat = trace.latitude;
	data(kk).slon = trace.longitude;
	data(kk).stnm = trace.station;
    data(kk).elev = trace.elevation;
	data(kk).timeaxis = new_timeaxis;
	data(kk).dataZ = dataZ;
	data(kk).dataR = dataR;
	data(kk).dataT = dataT;
    data(kk).dataZfilt = filtfilt(fb,fa,dataZ);
    data(kk).dataRfilt = filtfilt(fb,fa,dataR);
    data(kk).dataTfilt = filtfilt(fb,fa,dataT);
end
ntraces = length(data);

% pick component to look at
[data.dat] = data.(['data' component]);
[data.datfilt] = data.(['data' component 'filt']);

% pick low snr
noisy = [];
tt = data(1).timeaxis;
js = find(tt >= xcorwind(1) & tt < xcorwind(2));
jn = [2*js(1)-js(end)-1:js(1)-1]';
for is = 1:ntraces
    dat = data(is).datfilt;
    ds = detrend(dat(js),'constant');
    dn = detrend(dat(jn),'constant');
    snr = var(ds)/var(dn);
    if snr < snrmin, noisy = [noisy;is]; end
end
data(noisy) = []; % delete low snr data
ntraces = length(data);

%% cross correlate on travel time
% put into one matrix
tt = data(1).timeaxis;
inds = find(tt >= xcorwind(1) & tt < xcorwind(2));
for is = 1:ntraces
    zz(:,is) = data(is).datfilt(inds)./max(data(is).datfilt);
end

% cross correlate to get all aligned on arrival
% recursively xcor untill there are no traces w/ low acor  
acor = 0;
gd  = 1:ntraces;
% kill nans
gd(isnan(mean(zz))) = [];
while any(acor<acormin)
fprintf('Cross-correlating %.0f records\n',length(gd))
[dcor,dcstd,dvcstd,acor]=xcortimes(zz(:,gd),re_dt,-xcorwind(1),xcorlagmax,1);
gd(acor<acormin) = [];
figure(gcf) 
end
fprintf('From %.0f stas, keeping %.0f with acor > %.2f\n',ntraces,length(gd),acormin) 
data = data(gd);
ntraces = length(data);

% shift traces to align
old_timeaxis = tt;
new_timeaxis = tt(1) + xcorlagmax : re_dt : tt(end) - xcorlagmax; 
for is = 1:ntraces
    data(is).timeaxis = new_timeaxis;
    data(is).dataZ = interp1(old_timeaxis-dcor(is),data(is).dataZ,new_timeaxis);
    data(is).dataR = interp1(old_timeaxis-dcor(is),data(is).dataR,new_timeaxis);
    data(is).dataT = interp1(old_timeaxis-dcor(is),data(is).dataT,new_timeaxis);
    data(is).dataZfilt = interp1(old_timeaxis-dcor(is),data(is).dataZfilt,new_timeaxis);
    data(is).dataRfilt = interp1(old_timeaxis-dcor(is),data(is).dataRfilt,new_timeaxis);
    data(is).dataTfilt = interp1(old_timeaxis-dcor(is),data(is).dataTfilt,new_timeaxis);
    data(is).dcor = dcor(is);
    data(is).acor = acor(is);
end

% pick component to look at
[data.dat] = data.(['data' component]);
[data.datfilt] = data.(['data' component 'filt']);

stack = zeros(size(data(1).datfilt));
figure(2), clf, hold on
for is = 1:length(data)
%     plot(data(is).timeaxis - dcor(is),data(is).odataZfilt)
    hp = plot(data(is).timeaxis ,data(is).datfilt);
    stack = stack+data(is).datfilt;
    set(hp,'color',colour_get(data(is).slon,max([data.slon]),min([data.slon])))
end
xlim([-50 50])
stack = stack/ntraces;
plot(data(1).timeaxis,stack,'r','LineWidth',2)

%% Plot stations to be used - only a few
% only use data from swath of stations
gd = [data.slat] < latlim(2) & [data.slat] > latlim(1) & [data.slon] < lonlim(2) & [data.slon] > lonlim(1);
data = data(gd);
% eval(sprintf('for is = 1:length(data), data(is).dat = data(is).data%s; end',component))

%% work out distance to ridge
roughjdf = [-130.5,44.36;-128.7,48.98];
Xrdg = dist2line(roughjdf(1,:),roughjdf(2,:),[[data.slon]',[data.slat]'])'; 
% [data.Xrdg] = Xrdg << not sure why this doesn't work
for is = 1:length(data)
    data(is).Xrdg = Xrdg(is);
end

figure(3), clf, hold on
set(gcf,'pos',[ 10  50  900  680])
CascMAP
scatter([data.slon],[data.slat],70,colour_get(abs(Xrdg),max(abs(Xrdg)),min(abs(Xrdg))),'filled')
text([data.slon],[data.slat]+0.1,{data.stnm})

% Reference station
irst = find(strcmp(refsta,{data.stnm}));
plot(data(irst).slon, data(irst).slat,'ok','MarkerSize',10,'LineWidth',2.5)


%% calculate spectra

samprate = 1./unique(round_level(diff(data(1).timeaxis),0.001));
nstas = length(data);
dt = 1./samprate;
nyq=0.5./dt;

tt = data(1).timeaxis;
jn = find(tt >= t0-wlen & tt < t0);
js = find(tt >= t0 & tt < t0 + wlen);
nw=length(jn);
wdo=tukeywin(nw,0.2)';   % 0.2 gives 10% cosine taper


figure(4), clf
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
figure(5), clf, hold on
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
 
%% plot spectral ratios and measure tstar
figure(6), clf, hold on
for is = 1:nstas
if ~any(is==gd), continue, end    
ind = frq < hifrq;
hp = plot(frq(ind),data(is).lnR(ind),'o-');
set(hp,'color',colour_get(abs(Xrdg(is)),max(abs(Xrdg)),min(abs(Xrdg))))
xlim([0 hifrq])
fo = fit(frq(ind),data(is).lnR(ind),'poly1');
data(is).dtstar = fo.p1/(-pi);
data(is).dtstar_cb = confint(fo); data(is).dtstar_cb = data(is).dtstar_cb(:,1);
end


%% plot delta tstar and delta t
figure(7), clf, set(gcf,'pos',[59   258   871   513])
% delta tstar
subplot(2,1,1), hold on
plot(Xrdg(gd),[data(gd).dtstar],'.','MarkerSize',30)
text(Xrdg(gd),[data(gd).dtstar]+0.02,{data(gd).stnm})
% eb = [data(gd).dtstar_cb]';
% errorbar([data(gd).slon]',[data(gd).dtstar]',eb(:,1),eb(:,2),'o')
grid on
ylabel('$\Delta t^{\ast}$ (+ive means lower Q)','interpreter','latex','Fontsize',15)

% delta t
subplot(2,1,2), hold on
plot(Xrdg(gd),[data(gd).dcor],'.','MarkerSize',30)
grid on
ylabel('$\Delta t$ (+ive means later/slower)','interpreter','latex','Fontsize',15)



figure(4)
for is = 1:nstas
if ~any(is==gd), continue, end    
%% spectral plot
subplot(211), hold on
n4=length(find(frq<=0.8*nyq));
hps = plot(frq(1:n4),log10(data(is).specss(1:n4)),'b','LineWidth',1.5); % << here using smoothed spectrum 
hpn = plot(frq(1:n4),log10(data(is).specn(1:n4)),'k','LineWidth',1.5);
set(hps,'color',colour_get(abs(Xrdg(is)),max(abs(Xrdg)),min(abs(Xrdg))))
xlabel('Hz'); 
ylabel('amplitude, nm/Hz')
xlim([0 hifrq]);
ylim([-10,-3]);
end

end % if redo

%% find stas within mindist of section

isgdobs = [data(gd).elev]<0;
kpd = 78; % km per degree of longitude


%% Final proposal plot with everything on same axes
%% plot delta tstar and delta t
figure(17), clf, set(gcf,'pos',[59   258   871   613])
%% topo
subplot(5,1,1), hold on
plot(section_x,section_z,'k','LineWidth',2)
% figure things
xlim([-80 720]),ylim([-5001 5001])
set(gca,'visible','off')
% draw axes
plot([-80 720],[0 0],'--k') % x-axis
plot([-80 -80],[-5000 5000],'k','LineWidth',2)
plot([-80 -72],[-5000 -5000],'k',[-80 -72],[5000 5000],'k','LineWidth',1)
text(-80,7000,'5000','Fontsize',10,'HorizontalAlignment','center','VerticalAlignment','middle')
text(-80,-7000,'-5000','Fontsize',10,'HorizontalAlignment','center','VerticalAlignment','middle')
text(-150,0,'Elev (m)','Fontsize',13,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','Rotation',90)
% title
text(300,-7000,sprintf('%s $~$  Comp: %s, $f_{hi}~$: %.1f',datestr(evtime,'yyyy-mm-dd HH:MM:SS'),component,hifrq),'Fontsize',14,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex')
% ridge axis label
plot(0,4500,'vk','MarkerSize',6,'MarkerFaceColor','k')
text(0,6500,'Axis','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'Interpreter','latex')
% deformation front label
plot(305,4500,'vk','MarkerSize',6,'MarkerFaceColor','k')
text(305,6500,'DF','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'Interpreter','latex')
% Coastline label
plot(440,4500,'vk','MarkerSize',6,'MarkerFaceColor','k')
text(440,6500,'Coastline','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'Interpreter','latex')
% Arc label
plot(616,4500,'vk','MarkerSize',6,'MarkerFaceColor','k')
text(616,6500,'Arc','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'Interpreter','latex')


%% delta tstar
subplot(5,1,2:3), hold on
plot(Xrdg(gd( isgdobs))*kpd,[data(gd (isgdobs)).dtstar],'.b','MarkerSize',18)
% plot(Xrdg(gd(~isgdobs))*kpd,[data(gd(~isgdobs)).dtstar],'.r','MarkerSize',18)
plot(([data(gd(~isgdobs)).slon]+124)*kpd + 440,[data(gd(~isgdobs)).dtstar],'.r','MarkerSize',18)
% figure things
grid on
text(-150,0.5,'$\Delta t^{\ast}_S$','Fontsize',13,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','Rotation',90)
xlim([-80 720])

%% delta t
subplot(5,1,4:5), hold on
plot(Xrdg(gd( isgdobs))*kpd,[data(gd (isgdobs)).dcor],'.b','MarkerSize',18)
% plot(Xrdg(gd(~isgdobs))*kpd,[data(gd(~isgdobs)).dcor],'.r','MarkerSize',18)
plot(([data(gd(~isgdobs)).slon]+124)*kpd + 440,[data(gd(~isgdobs)).dcor],'.r','MarkerSize',18)
% figure things
grid on
text(-150,0,'$\Delta t_S$','Fontsize',13,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','Rotation',90)
text(320,-4.5,'Distance from ridge (km)','Fontsize',13,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex')
xlim([-80 720]), ylim([-3 3])

save2pdf(17,'observed_dt_dtstar','figs');

%% Final proposal map with everything on same axes
figure(3), clf, hold on
set(gcf,'pos',[ 10  50  900  680])
CascMAP
% scatter([data(gd).slon],[data(gd).slat],70,colour_get(abs(Xrdg(gd)),max(abs(Xrdg)),min(abs(Xrdg))),'filled')
scatter([data(gd( isgdobs)).slon],[data(gd( isgdobs)).slat],70,'b','filled')
scatter([data(gd(~isgdobs)).slon],[data(gd(~isgdobs)).slat],70,'r','filled')
text([data(gd).slon],[data(gd).slat]+0.15,{data(gd).stnm},'FontSize',7)
box on

% Reference station
irst = find(strcmp(refsta,{data.stnm}));
plot(data(irst).slon, data(irst).slat,'ok','MarkerSize',10,'LineWidth',2.5)

save2png(3,'proposal_mapview','figs');

%% Final proposal waveforms with everything on same axes
Proposal_fig3
