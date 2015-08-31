% script to go through data and reset to raw traces and channels=
clear all
addpath('matguts')

overwrite = true;

% antelope db details
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascattendb';

% path to top level of directory tree for data
datadir = '~/Work/CASCADIA/DATA/'; % needs final slash

%% get to work
db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,evtimes] = dbgetv(dbor,'orid','time');
norids = dbnrecs(dbor);
dbclose(db);

for ie = 85:85 % 1:norids % loop on orids
    fprintf('\n Orid %.0f %s \n\n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
    evdir = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
    load([datadir,evdir,'_datinfo.mat'])
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end

    for is = 1:length(datinfo) % loop on stas

        sta = datinfo(is).sta; % sta name
        load([datadir,evdir,sta,'.mat']); % load sta data for this evt
        fprintf('Station %-5s...',sta)
        data.dat    = data.raw.dat;
        data.chans  = data.raw.chans;
        data.NEZ    = false;
        data.rmresp = false;
        data.rmtilt = false;
        data.rmcomp = false;
        
        [data.phases.xtime]   = deal([]);
        [data.phases.xartime] = deal([]);
        [data.phases.xacor]   = deal([]);
        
        datinfo(is).chans  = data.raw.chans.component;
        datinfo(is).NEZ    = false;
        datinfo(is).rmresp = false;
        datinfo(is).rmtilt = false;
        datinfo(is).rmcomp = false;
        
        %% log reset and save
        fprintf(' reset\n')
        save([datadir,evdir,sta],'data')
        save([datadir,evdir,'_datinfo'],'datinfo')
            
    end % loop on stas

	%% sum up
    fprintf(' STA  CHAN  NEZ  resp  tilt  comp\n')
    for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f     %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).rmtilt,datinfo(is).rmcomp); end

end
