% script to build body wave dataset
% uses antelope origin table and site(chan) tables, looping through and
% then using the IRISrequest tools to get and store the data in a directory
% tree where each event has a folder containing .mat files that are the
% data for each staition.
% clear all
close all
cd ~/Documents/MATLAB/CASC_atten
% mount_drive('DATA','zeilon','eilon.ldeo.columbia.edu')

datawind = [-100 1700]; % time window in seconds after event to [start end]

phases = 'P,S,PKS,SKS';

overwrite = false;

resamprate = 5; % leave empty or zero to use existing samprate

getnoise = true;
OBSnoiseprewind = [-43200 0]; % time window in seconds after event to [start end] <== 12 hours in advance

% antelope db details
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';

% path to top level of directory tree for data
datadir = '/Volumes/DATA_mini/CASCADIA/DATA/'; % needs final slash
% datadir = '~/Work/CASCADIA/DATA/';

javaaddpath('/Users/zeilon/Documents/MATLAB/IRIS-WS-2.0.15.jar')
javaaddpath('IRIS-WS-2.0.15.jar')
addpath('matguts')


%% get event details
db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,elats,elons,edeps,evtimes] = dbgetv(dbor,'orid','lat','lon','depth','time');
norids = dbnrecs(dbor);
dbclose(db);

%% get station details
db = dbopen([dbdir,dbnam],'r');
dbsi = dblookup_table(db,'site');
[stas,slats,slons,selevs,nwk,statype] = dbgetv(dbsi,'sta','lat','lon','elev','refsta','statype');
nstas = dbnrecs(dbsi);
dbclose(db);

for ie = 220:220 % 1:norids
    % sort out event stuff
    orid = orids(ie);
    elat = elats(ie); elon = elons(ie); edep = edeps(ie); 
    evtime = evtimes(ie);
    evdir = [num2str(orid,'%03d'),'_',epoch2str(evtime,'%Y%m%d%H%M')];
    
    % make event directory
    if exist([datadir,evdir],'dir')~=7, mkdir(datadir,evdir); end
    
    % datinfo
    datinfofile = [datadir,evdir,'/_datinfo'];
    if exist([datinfofile,'.mat'],'file')~=2 || overwrite==true
        datinfo = struct('sta',[],'chans',[],'NEZ',false,'rmresp',false,'rmtilt',false,'rmcomp',false,'spectra',false);
    else
        load(datinfofile);
    end
        
    
    fprintf('REQUESTING DATA FOR EVENT %.0f (%s)\n',orid,evdir)
    for is = 259:nstas % 1:nstas
        sta = stas{is};
        datafile = [datadir,evdir,'/',stas{is}];
        
        if overwrite==false && exist([datafile,'.mat'],'file')==2, continue; end
        
        % this is where we'll pull out the channel and station on/off info
        db = dbopen([dbdir,dbnam],'r');
        dbsch = dblookup_table(db,'sitechan');
        dbschs = dbsubset(dbsch,sprintf('sta == "%s"',stas{is}));
        nchans = dbnrecs(dbschs);
        if nchans==0, dbclose(db); continue; end
        [chans,ondates,offdates,azimuths,dips] = dbgetv(dbschs,'chan','ondate','offdate','hang','vang');
        dbclose(db);
        if ~iscell(chans), chans = {chans}; end;
        
        % calc. data window
        if strcmp(statype(is),'OBS') % if OBS, if the option is selected, grab big noise window too!
            if getnoise
                STARTtime = evtime + datawind(1) + OBSnoiseprewind(1) - 1; % inc buffer
                ENDtime   = evtime + datawind(2) + OBSnoiseprewind(2) + 1; % inc buffer
            else
                STARTtime = evtime + datawind(1) - 1; % inc buffer
                ENDtime   = evtime + datawind(2) + 1; % inc buffer
            end
        else % if LAND
        STARTtime = evtime + datawind(1) - 1; % inc buffer
        ENDtime   = evtime + datawind(2) + 1; % inc buffer
        end
        
        STARTtime_str = epoch2str(STARTtime,'%Y-%m-%d %H:%M:%S.%s');
        ENDtime_str   = epoch2str(ENDtime,'%Y-%m-%d %H:%M:%S.%s');
        
        
        % check this station was alive for this event
        ondate = num2str(unique(ondates)); offdate = num2str(unique(offdates));
        onepoch = str2epoch([ondate(1:4),'/',ondate(5:7),' 00:00:00']);
        offepoch=str2epoch([offdate(1:4),'/',offdate(5:7),' 23:59:59']);
        if onepoch  > (STARTtime) , continue, end % continue if stat turned on after datwind start
        if offepoch < (ENDtime), continue, end % continue if stat turned off before datwind end
        
        % make string list of chans to request
        chreq = []; for ic = 1:length(chans), chreq = [chreq,chans{ic},',']; end; chreq = chreq(1:end-1); %#ok<AGROW>
        
        % exceptional names
        if strcmp(sta,'M02CO') || strcmp(sta,'M04CO') || strcmp(sta,'M02CL') || strcmp(sta,'M04CL')
            sta_use = sta(1:4);
        else
            sta_use = sta;
        end
        fprintf('   request station %.0f %s... ',is,stas{is})
        
        % ======== GET THE DATA =======
        trace=irisFetch.Traces(nwk{is},sta_use,'*',chreq,STARTtime_str,ENDtime_str);
        
        if isempty(trace), fprintf('NO DATA\n'); continue; end
        [ trace ] = fixtrace( trace );
        fprintf('got chans '); for ic = 1:length(trace), fprintf('%s, ',trace(ic).channel); end
        % if B and H channels - only keep H (might need higher samprate?)
        if strcmp([trace.channel],'BHEBHNBHZHHEHHNHHZ'), trace = trace(4:6); end
        if strcmp([trace.channel],'HHEHHNHHZBHEBHNBHZ'), trace = trace(1:3); end
        
        % station details from trace
        sta_dts = struct('name',trace(1).station,...
                         'slat',trace(1).latitude,...
                         'slon',trace(1).longitude,...
                         'selev',trace(1).elevation);
        
        % channels' details
        comp = cell(1,length(trace)); for ic = 1:length(trace), comp{ic} = trace(ic).channel(end); end
        chan_dts = struct('name',{{trace.channel}},...
                          'component',{comp},...
                          'sensitivity',[trace.sensitivity],...
                          'dip',[trace.dip],...
                          'azimuth',[trace.azimuth]);
                      
        [gcarc,esaz] = distance(elat,elon,sta_dts.slat,sta_dts.slon);
        [~,seaz] = distance(sta_dts.slat,sta_dts.slon,elat,elon);   

        % SAMPRATE
        if resamprate
            samprate = resamprate;
        else
            samprate = round(unique([trace.sampleRate]));
            if length(samprate)>1, 
                fprintf(' different samprates, downsamp to min'); 
                samprate = round(min(unique([trace.sampleRate])));
            end
        end
        
      
        % TIME
        tt0 = STARTtime + 1;
        tt1 = ENDtime - 1 - 1./samprate;
        tt = [tt0:(1./samprate):tt1]'; 
        
        % SAFETY
        if any([trace.startTime]>epoch2serial(tt0))
            fprintf('REQUESTED DATA ONLY STARTS AFTER DESIRED WINDOW START!!\n')
            continue            
        end
        if any([trace.endTime]<epoch2serial(tt1))
            fprintf('REQUESTED DATA ENDS BEFORE DESIRED WINDOW END!!\n')
            continue
        end
        
        % data
        nsamps = length(tt);
        dat = zeros(nsamps,length(trace));
        for id = 1:length(trace)
        dat(:,id) = interp1(linspace(serial2epoch(trace(id).startTime),...
                                     serial2epoch(trace(id).endTime),...
                                     trace(id).sampleCount),...
                                     trace(id).data,    tt);
        end
               
        % fix nans
        nandat = find(isnan(dat));
        if ~isempty(nandat)
            if length(nandat) > 2*length(trace)
                fprintf('lots of nans - look out')
            end
            fprintf('fixing nans')
            dat(nandat) = 0;
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
                       'phases',TT, 'dat',dat,'tt',tt,...
                       'raw',struct('chans',chan_dts,'dat',dat),...
                       'NEZ',false,'rmresp',false,'rmtilt',false,'rmcomp',false);

        datinfo(is,1) = struct('sta',stas{is},'chans',{chan_dts.component},'NEZ',false,'rmresp',false,'rmtilt',false,'rmcomp',false,'spectra',false);
        % save
        save(datafile,'data')
        save(datinfofile,'datinfo')
        fprintf('\n')
    end %loop on stas
    kill = []; for ii = 1:length(datinfo), if isempty(datinfo(ii).sta), kill = [kill;ii]; end; end
    datinfo(kill) = [];
    save(datinfofile,'datinfo')
end % loop on evts


