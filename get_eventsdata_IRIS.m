% script to build body wave dataset
% uses antelope origin table and site(chan) tables, looping through and
% then using the IRISrequest tools to get and store the data in a directory
% tree where each event has a folder containing .mat files that are the
% data for each staition.

datawind = [-200 1800]; % time window in seconds after event to [start end]
phases = 'P,S,PKS,SKS';
overwrite = false;

% antelope db details
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascattendb';

% path to top level of directory tree for data
datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash

javaaddpath('/Users/zeilon/Documents/MATLAB/IRIS-WS-2.0.14.jar')


%% get event details
db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,elats,elons,edeps,evtimes] = dbgetv(dbor,'orid','lat','lon','depth','time');
norids = dbnrecs(dbor);
dbclose(db);

%% get station details
db = dbopen([dbdir,dbnam],'r');
dbsi = dblookup_table(db,'site');
[stas,slats,slons,selevs,nwk] = dbgetv(dbsi,'sta','lat','lon','elev','refsta');
nstas = dbnrecs(dbsi);
dbclose(db);

for ie = 61:85 % 1:norids
    % sort out event stuff
    orid = orids(ie);
    elat = elats(ie); elon = elons(ie); edep = edeps(ie); 
    evtime = evtimes(ie);
    evdir = [num2str(orid,'%03d'),'_',epoch2str(evtime,'%Y%m%d%H%M')];
    
    % make event directory
    if exist([datadir,evdir],'dir')~=7, mkdir(datadir,evdir); end
    
    % calc. data window
    waveform_start_time = epoch2str(evtime + datawind(1) - 1,'%Y-%m-%d %H:%M:%S'); % inc buffer
    waveform_end_time   = epoch2str(evtime + datawind(2) + 1,'%Y-%m-%d %H:%M:%S'); % inc buffer
    
    fprintf('REQUESTING DATA FOR EVENT %.0f (%s)\n',orid,evdir)
    for is = 1:nstas
        sta = stas{is};
        if overwrite==false && exist([datadir,evdir,'/',stas{is},'.mat'],'file')==2, continue; end
        
        % this is where we'll pull out the channel and station on/off info
        db = dbopen([dbdir,dbnam],'r');
        dbsch = dblookup_table(db,'sitechan');
        dbschs = dbsubset(dbsch,sprintf('sta == "%s"',stas{is}));
        nchans = dbnrecs(dbschs);
        if nchans==0, dbclose(db); continue; end
        [chans,ondates,offdates,azimuths,dips] = dbgetv(dbschs,'chan','ondate','offdate','hang','vang');
        dbclose(db);
        if ~iscell(chans), chans = {chans}; end;
        
        % check this station was alive for this event
        ondate = num2str(unique(ondates)); offdate = num2str(unique(offdates));
        onepoch = str2epoch([ondate(1:4),'/',ondate(5:7),' 00:00:00']);
        offepoch=str2epoch([offdate(1:4),'/',offdate(5:7),' 23:59:59']);
        if onepoch  > evtime + datawind(1), continue, end % continue if stat turned on after datwind start
        if offepoch < evtime - datawind(2), continue, end % continue if stat turned off before datwind end
        
        % make string list of chans to request
        chreq = []; for ic = 1:length(chans), chreq = [chreq,chans{ic},',']; end; chreq = chreq(1:end-1); %#ok<AGROW>
        
        % exceptional names
        if strcmp(sta,'M02CO') || strcmp(sta,'M04CO') || strcmp(sta,'M02CL') || strcmp(sta,'M04CL')
            sta = sta(1:4);
        end
        
        fprintf('   request station %.0f %s... ',is,stas{is})
        trace=irisFetch.Traces(nwk{is},sta,'*',chreq,waveform_start_time,waveform_end_time);
        if isempty(trace), fprintf('NO DATA\n'); continue; end
        [ trace ] = parse_trace( trace );
        fprintf('got chans '); for ic = 1:length(trace), fprintf('%s, ',trace(ic).channel); end, fprintf('\n')
        
        % if B and H channels - only keep H (might need higher samprate?)
        if strcmp([trace.channel],'BHEBHNBHZHHEHHNHHZ'), trace = trace(1:3); end
        if strcmp([trace.channel],'HHEHHNHHZBHEBHNBHZ'), trace = trace(4:6); end
        
        % station details from trace
        sta_dts = struct('name',trace(1).station,...
                         'slat',trace(1).latitude,...
                         'slon',trace(1).longitude,...
                         'selev',trace(1).elevation);
        
        % channels' details
        for ic = 1:length(trace), comp{ic} = trace(ic).channel(end); end
        chan_dts = struct('name',{{trace.channel}},...
                          'component',comp,...
                          'sensitivity',[trace.sensitivity],...
                          'dip',[trace.dip],...
                          'azimuth',[trace.azimuth],'rot',false);
                      
        [gcarc,esaz] = distance(elat,elon,sta_dts.slat,sta_dts.slon);
        [~,seaz] = distance(sta_dts.slat,sta_dts.slon,elat,elon);   
        samprate = round(unique([trace.sampleRate])); 
        if length(samprate)>1, error('differnt samprates\n'); end
      
        % time
        startTime = evtime + datawind(1);
        endTime   = evtime + datawind(2) - 1./samprate;
        tt = [startTime:(1./samprate):endTime]'; 
        
        % safety
        if any([trace.startTime]>epoch2serial(startTime))
            error('Requested data only starts after desired window start')
        end
        if any([trace.endTime]<epoch2serial(endTime))
            error('Requested data ends before desired window end')
        end
        
        % data
        nsamps = abs(diff(datawind))*samprate;
        dat = nan(nsamps,length(trace));
        for id = 1:length(trace)
        dat(:,id) = interp1(linspace(serial2epoch(trace(id).startTime),...
                                     serial2epoch(trace(id).endTime),...
                                     trace(id).sampleCount),...
                                     ...
                                     trace(id).data,    tt);
        end
                
        % phase and distance details
        TT = tauptime('event',[elat,elon],'depth',edep,'station',[sta_dts.slat,sta_dts.slon],'phases',phases );
        for ip = 1:length(TT)
            TT(ip).artime = TT(ip).time + evtime;
        end
        
        % finally put into structure
        data = struct('station',sta_dts,'network',trace(1).network,...
                       'chans',chan_dts,'samprate',samprate,'nsamps',nsamps,...
                       'gcarc',gcarc,'seaz',seaz,'esaz',esaz,...
                       'phases',TT, 'dat',dat,'tt',tt);

        % save
        save([datadir,evdir,'/',stas{is}],'data')
        
    end
end


