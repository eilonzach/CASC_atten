clear all
%% find event
% gdevts_fn([2013,08,08],[2014,05,14],[6.2 7.1],[],[],[],[46 -125.5],[20 90]);
% return

javaaddpath('IRIS-WS-2.0.6.jar')

%% Data inputs
time_window = [200 200]; %time (s) [before after] predicted phase to window
phase = 'S';
%% Station inputs
network = '*'; %'*' for all, else choose '7D' OBS, or WA, TA, 
chans = {'BH?','HH?'};
latlim = [43 49];
lonlim = [-131 -120];



%% Event inputs
%evtime format: 'yyyy-mm-dd HH:MM:SS'
% evtime = '2014-01-01 16:03:29'; magnitude = 6.5; % event 1
% evtime = '2014-06-23 20:53:09'; magnitude = 7.9; % event 2
% evtime = '2013-10-25 17:10:20'; magnitude = 7.1; % event 3 - Japan
% evtime = '2013-09-25 16:42:20'; magnitude = 7.1; % event 4 - Peru/Chile
% evtime = '2014-05-13 06:35:24'; magnitude = 6.5; % event 5 - Nicaragua
% evtime  = '2013-09-04 00:18:23'; magnitude = 6.5; % event 6 - SOUTH OF HONSHU, JAPAN dep=402.0  delta=72.6   baz=295  
evtime  = '2013-10-01 03:38:22'; magnitude = 6.7; % event 7 - SEA OF OKHOTSK dep=573.0  delta=50.5  baz=310  


search_window = 1; %time before and after evtime to search, in hrs


%% grab some data
load('coast')
jdf = dlmread('mapdata/ridge_xy');
fzs = dlmread('mapdata/transforms_xy');

%% Get event details
search_wstart = datenum(evtime,'yyyy-mm-dd HH:MM:SS') - search_window/24;
search_wend   = datenum(evtime,'yyyy-mm-dd HH:MM:SS') + search_window/24;
search_wstart = datestr(search_wstart,'yyyy-mm-dd HH:MM:SS');
search_wend = datestr(search_wend,'yyyy-mm-dd HH:MM:SS');

evinfo = irisFetch.Events('boxcoordinates',[-90 90 -180 180],...
		'startTime',search_wstart,'endTime',search_wend,...
		'MinimumMagnitude',magnitude-0.1,'maximumMagnitude',magnitude+0.1);
if length(evinfo)>1, error('Found multiple events, refine search'), end

%% Get station details
evtime = evinfo.PreferredTime;
start_time = evtime;
end_time = datenum(start_time,'yyyy-mm-dd HH:MM:SS') + 2/24; % in days
end_time = datestr(end_time,'yyyy-mm-dd HH:MM:SS');

stinfo = struct([]);
for ic = 1:length(chans)
junk = irisFetch.Stations('channel',network,'*','*',chans{ic},...
					'startTime',start_time,'endTime',end_time,...
                    'boxcoordinates',[latlim, lonlim]);
stinfo = [stinfo; junk'];
end

nstas = length(stinfo);
fprintf('Downloading data for %.0f stations\n',nstas)


% plot stations trying to get data from
figure(9), clf, hold on
scatter([stinfo.Longitude],[stinfo.Latitude],50,...
    colour_get([stinfo.Longitude],max([stinfo.Longitude]),min([stinfo.Longitude])),'filled')
plot(long,lat)
plot(jdf(:,1),jdf(:,2),'k','Linewidth',2)
plot(fzs(:,1),fzs(:,2),'--k','Linewidth',1.5)
xlim(lonlim)
ylim(latlim)
title('Stations requesting data from','FontSize',14)
pause(0.01)

%% Get data, parse into structure
traces = [];
notr = [];
for is = 1:nstas
    network = stinfo(is).NetworkCode;
    sta_nam = stinfo(is).StationCode;
    TT = tauptime('event',[evinfo.PreferredLatitude,evinfo.PreferredLongitude],...
                  'depth',evinfo.PreferredDepth,...
                  'station',[stinfo(is).Latitude,stinfo(is).Longitude],...
                  'phases',phase );
    stinfo(is).Ptime = TT.time;
    stinfo(is).twin  = time_window;
    t1 = TT.time-time_window(1)-5; % time from evtime in s to start window (w 5s buffer)
    t2 = TT.time+time_window(2)+5; % time from evtime in s to start window (w 5s buffer)
    waveform_start_time = datestr(datenum(evtime,'yyyy-mm-dd HH:MM:SS') + t1/86400,'yyyy-mm-dd HH:MM:SS');
    waveform_end_time   = datestr(datenum(evtime,'yyyy-mm-dd HH:MM:SS') + t2/86400,'yyyy-mm-dd HH:MM:SS');

    fprintf('Downloading traces for %s from %s to %s (%.0f/%.0f)\n',sta_nam,waveform_start_time,waveform_end_time,is,nstas)

    try
    trace=irisFetch.Traces(network,sta_nam,'*','BH?',...
                                waveform_start_time,waveform_end_time,'includePZ');
    catch e
        
    try 
    trace=irisFetch.Traces(network,sta_nam,'*','HH?',...
                                waveform_start_time,waveform_end_time,'includePZ');
    catch e
        notr = [notr,is]; continue
    end
    end
    
    if isempty(trace),  notr = [notr,is]; continue; end
    
%     if length(trace) > 3
%         if strcmp(trace(1).channel,trace(4).channel) & ~strcmp(trace(1).location,trace(4).location)
%             if trace(1).sampleRate > trace(4).sampleRate
%                 trace = trace(1:3);
%             else
%                 trace = trace(4:6);
%             end
%         else
%         error('more to trace than thought')
%         junk = trace(1); 
%         junk.sampleCount = sum([trace.sampleCount]);
%         junk.data = [junk.data;trace(2).data];
%         junk.endTime = trace(2).endTime;
%         trace = junk;
%         end
%     end
    trace = fixtrace(trace);
    
    for ic = 1:3
    cha_nam = trace(ic).channel(1:3);   
    %% get orientations:
    % ideally do this through antelope but license expired
    fid = fopen('cascyr3.sitechan','r');
    A = textscan(fid,'%6s %6s   %7.0f      %3.0f  %7.0f %s %6.4f %5.1f %4.1f %s %f');
    fclose(fid);
    newaz = A{8}(strcmp(A{1},sta_nam) & strcmp(A{2},cha_nam));
    newdip= A{9}(strcmp(A{1},sta_nam) & strcmp(A{2},cha_nam)) - 90;
    if  ~isempty(newaz)
        trace(ic).azimuth = newaz;
        trace(ic).dip     = newdip;
    end
    
    end % loop on chans
    
	if isempty(traces), traces = trace; else traces(is,:) = trace; end
    
end

    stinfo(notr) = [];
    traces(notr(notr<=length(traces)),:) = [];
    
% plot stations trying to get data from
figure(10), clf, hold on
scatter([stinfo.Longitude],[stinfo.Latitude],50,...
    colour_get([stinfo.Longitude],max([stinfo.Longitude]),min([stinfo.Longitude])),'filled')
plot(long,lat)
plot(jdf(:,1),jdf(:,2),'k','Linewidth',2)
plot(fzs(:,1),fzs(:,2),'--k','Linewidth',1.5)
xlim(lonlim)
ylim(latlim)
title('Stations got data from','FontSize',14)
pause(0.01)

%% Save for posterity
save('eg_evt_sta_info.mat','evinfo','stinfo')
save('eg_traces.mat','traces')

evstr = sprintf('_%s_%s',datestr(datenum(evtime,'yyyy-mm-dd HH:MM:SS'),'yyyymmdd'),phase);
copyfile('eg_evt_sta_info.mat',['eg_evt_sta_info' evstr '.mat']);
copyfile('eg_traces.mat',['eg_traces' evstr '.mat']);

% Fetch data form IRIS
%   args: Net, Sta, Loc, Cha, Starttime, Endtime [,quality][,includePZ][,verbosity]
