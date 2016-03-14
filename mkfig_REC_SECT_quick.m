% quick and dirty script to make a record section for a given event
clear all
orid = 188;
comp = 'Z';

amp = 0.1;
filtfs = [0.04 2]; % in hz

% antelope db details
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';

% path to top level of directory tree for data
datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash
% datadir = '~/Work/CASCADIA/DATA/';


%% get to work
db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
dbors = dbsubset(dbor,sprintf('orid == %.0f',orid));
[orids,elats,elons,edeps,evtimes,mag] = dbgetv(dbors,'orid','lat','lon','depth','time','ms');
dbclose(db);
evdir = [num2str(orids,'%03d'),'_',epoch2str(evtimes,'%Y%m%d%H%M'),'/'];
datinfofile = [datadir,evdir,'_datinfo.mat'];
    
% prep file download info
fprintf('Looking in %s \n',[datadir,evdir])
load(datinfofile);
% if isempty(stafiles), error('No station mat files for this event'); end
nstas = length(datinfo);
figure(1), clf, set(gcf,'position',[10 10 1400 1000])
hold on

for is = 1:nstas;%nstas
    load([datadir,evdir,datinfo(is).sta,'.mat'])
    if ~data.NEZ, continue, end

    if any(strcmp(comp,{'E','N','Z'}))
        ic = find(strcmp(data.chans.component,comp),1,'first'); if isempty(ic), continue; end
        dat = data.dat(:,ic);
    elseif  any(strcmp(comp,{'R','T'}))
%         datZ = data.dat(:,strcmp(data.chans.component,'Z'));
        datN = data.dat(:,strcmp(data.chans.component,'N'));
        datE = data.dat(:,strcmp(data.chans.component,'E'));

        foraz = mod(data.seaz+180,360);
        datR =  datN*sind(foraz) + datE*sind(foraz);
        datT = -datN*sind(foraz) + datE*cosd(foraz);
        if strcmp(comp,'R'), dat = datR; end
        if strcmp(comp,'T'), dat = datT; end
    end
    
    % bit of cleaning
    dat = detrend(dat);
    
    W = 2*filtfs./data.samprate;
    [fb, fa] = butter(2,W);
    dat = filtfilt(fb,fa,dat);
    
    dat = dat/abs(max(dat));
    
    text(data.tt(1)-evtimes-1,data.gcarc,data.station.name,'horizontalalignment','left')
    plot(data.tt-evtimes,amp*dat+data.gcarc)
    
%     if strcmp(data.station.name,'SAO'), plot(data.tt-evtimes,amp*dat+data.gcarc,'g','LineWidth',2), end
%     if strcmp(data.station.name,'BKS'), plot(data.tt-evtimes,amp*dat+data.gcarc,'g','LineWidth',2), end
%     if strcmp(data.station.name,'BOZ'), plot(data.tt-evtimes,amp*dat+data.gcarc,'g','LineWidth',2), end
%     if strcmp(data.station.name,'J35C'), plot(data.tt-evtimes,amp*dat+data.gcarc,'g','LineWidth',2), end
%     if strcmp(data.station.name,'J43C'), plot(data.tt-evtimes,amp*dat+data.gcarc,'g','LineWidth',2), end
%     if strcmp(data.station.name,'J44C'), plot(data.tt-evtimes,amp*dat+data.gcarc,'g','LineWidth',2), end
    if strcmp(data.station.name,'BOZ'), plot(data.tt-evtimes,amp*dat+data.gcarc,'g','LineWidth',2), end
    if strcmp(data.station.name,'BB330'), plot(data.tt-evtimes,amp*dat+data.gcarc,'g','LineWidth',2), end
    if strcmp(data.station.name,'BB330'), plot(data.tt-evtimes,amp*dat+data.gcarc,'g','LineWidth',2), end
    if strcmp(data.station.name,'BB130'), plot(data.tt-evtimes,amp*dat+data.gcarc,'g','LineWidth',2), end
    if strcmp(data.station.name,'BB350'), plot(data.tt-evtimes,amp*dat+data.gcarc,'g','LineWidth',2), end
    if strcmp(data.station.name,'BB480'), plot(data.tt-evtimes,amp*dat+data.gcarc,'g','LineWidth',2), end
    
    for ip = 1:length(data.phases), plot(data.phases(ip).time,data.gcarc,'.g','MarkerSize',20); end
end

lms = axis;
gcs = [lms(3):0.2:lms(4)]';
phs = {'P','PP','Pdiff','S','SS','Sdiff','PKS','SKS','PS'};
predar = nan(length(gcs),length(phs)); 
for ig = 1:length(gcs)
for ip = 1:length(phs)
TT = tauptime('depth',edeps,'phases',phs{ip},'deg',gcs(ig));
if isempty(TT), continue; end 
predar(ig,ip) = TT.time;
end
end

plot(predar,gcs*ones(size(phs)),'-r')
ind = round(length(gcs)/2);
for ip = 1:length(phs)
    text(predar(ind,ip),gcs(ind),phs{ip},'FontSize',25,'color','k')
end