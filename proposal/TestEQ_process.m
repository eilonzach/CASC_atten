clear all

% evtime = '2014-01-01 16:03:29'; magnitude = 6.5; % event 1
% evtime = '2014-06-23 20:53:09'; magnitude = 7.9; % event 2
% evtime = '2013-10-25 17:10:20'; magnitude = 7.1; % event 3 - Japan
% evtime = '2013-09-25 16:42:20'; magnitude = 7.1; % event 4 - Peru/Chile
% evtime = '2014-05-13 06:35:24'; magnitude = 6.5; % event 5 - Nicaragua
% evtime = '2013-09-04 00:18:23'; magnitude = 6.5; % event 6 - SOUTH OF HONSHU, JAPAN dep=402.0  delta=72.6   baz=295  
% evtime  = '2013-10-01 03:38:22'; magnitude = 6.7; % event 7 - SEA OF OKHOTSK dep=573.0  delta=50.5  baz=310  

event = '20140513';
phase = 'P';

rmresp = 1;

resamprate = 40 ; % new, common sample rate
component = 'Z'; %'Z', 'R', or 'T'
filtfs = [25 0.5];
xcorwind = [-10 30];
xcorlagmax = 3;
acormin = 0.6;
snrmin = 3;

%% ================================================================
%% ============== shouldn't need to alter below here ==============

addpath /Users/zeilon/Documents/MATLAB/jingle/recreadplot

load(sprintf('eg_evt_sta_info_%s_%s.mat',event,phase));
load(sprintf('eg_traces_%s_%s.mat',event,phase));

evtime = datenum(evinfo.PreferredTime,'yyyy-mm-dd HH:MM:SS');
re_dt = 1./resamprate;

kk = 0;
for is = 1:size(traces,1)
    if isempty(traces(is,1).network), continue, end
%     st = min([traces(is,1:3).startTime]);
%     et = max([traces(is,1:3).endTime]);
%     b = (st - evtime )*24*3600  - stinfo(is).Ptime;
% 	  e = (et - evtime )*24*3600  - stinfo(is).Ptime;
%     sr = unique([traces(is,1:3).sampleRate]);
%     dt = 1/sr;
%     sc = (e-b)*sr;
%     dat = zeros(sc,3);
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

% % % make sure all same length
% % nsamps = [];
% % for is = 1:ntraces, nsamps = min([nsamps,length(data(is).timeaxis)]); end
% % for is = 1:ntraces
% %     data(is).timeaxis = data(is).timeaxis(1:nsamps); 
% %     data(is).dataZ = data(is).dataZ(1:nsamps);
% %     data(is).dataR = data(is).dataR(1:nsamps);
% %     data(is).dataT = data(is).dataT(1:nsamps);
% %     data(is).dataZfilt = data(is).dataZfilt(1:nsamps);
% %     data(is).dataRfilt = data(is).dataRfilt(1:nsamps);
% %     data(is).dataTfilt = data(is).dataTfilt(1:nsamps);
% % end

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
%     figure(78),clf,plot(jn,dn,'r',js,ds,'b'); title(snr)
    if snr < snrmin, noisy = [noisy;is]; end
end
data(noisy) = []; % delete low snr data
ntraces = length(data);

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
pause 
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

% % pick component to look at
% eval(sprintf('for is = 1:length(data), data(is).dat = data(is).data%s; end',component))
% eval(sprintf('for is = 1:length(data), data(is).datfilt = data(is).data%sfilt; end',component))
% pick component to look at
[data.dat] = data.(['data' component]);
[data.datfilt] = data.(['data' component 'filt']);


stack = zeros(size(data(1).datfilt));
figure(1), clf, hold on
for is = 1:length(data)
%     plot(data(is).timeaxis - dcor(is),data(is).odataZfilt)
    hp = plot(data(is).timeaxis ,data(is).datfilt);
    stack = stack+data(is).datfilt;
    set(hp,'color',colour_get(data(is).slon,max([data.slon]),min([data.slon])))
end
xlim([-50 50])
stack = stack/ntraces;
plot(data(1).timeaxis,stack,'r','LineWidth',2)

save('eg_cleandata.mat','data')
save('eg_cleaninfo.mat','evinfo')

