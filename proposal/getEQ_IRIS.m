javaaddpath('IRIS-WS-2.0.9.jar')

ev = [];
stinfo

network = '*'; %'*' for all, else choose '7D' OBS, or WA, TA, 
chans = {'BH?','HH?'};
latlim = [43 49];
lonlim = [-131 -120];


for yr = 2011:1:2015
startTime = datenum(yr,1,1);
endTime = datenum(yr+1,1,1);

tmp = irisFetch.Events('MinimumMagnitude',6,...
    'radialcoordinates',[46 -125 135 30],...
    'starttime',datestr(startTime,'yyyy-mm-dd HH:MM:SS'),...
    'endtime',datestr(endTime,'yyyy-mm-dd HH:MM:SS'));

ev = [ev;tmp'];
end

nevts = length(ev);
ndo = 10;
evdos = unique(random('unid',nevts,ndo,1));

nstas = 0;

for ie = 1:length(evdos)
    evtime = ev(evdos(ie)).PreferredTime
    start_time = evtime;
    end_time = datenum(start_time,'yyyy-mm-dd HH:MM:SS') + 0.5/24; % in days
    end_time = datestr(end_time,'yyyy-mm-dd HH:MM:SS');

    for ic = 1:length(chans)
    stinfo = irisFetch.Stations('channel',network,'*','*',chans{ic},...
                        'startTime',start_time,'endTime',end_time,...
                        'boxcoordinates',[latlim, lonlim]);
    fprintf('%.0f stations with %s data\n',length(stinfo),char(chans{ic}))
    nstas = nstas+length(stinfo);
    end
    % 	trace=irisFetch.Traces(network,sta_nam,'*','BH?',...
%                                 waveform_start_time,waveform_end_time,'includePZ') 
end

fprintf('Mean # stas per evt is %.1f\n',nstas/ndo)
