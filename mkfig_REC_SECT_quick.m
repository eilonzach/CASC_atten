% quick and dirty script to make a record section for a given event
clear all
orid = 86;
comp = 'Z';

% antelope db details
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascattendb';

% path to top level of directory tree for data
datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash
datadir = '~/Work/CASCADIA/DATA/';


%% get to work
db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
dbors = dbsubset(dbor,sprintf('orid == %.0f',orid));
[orids,elats,elons,edeps,evtimes,mag] = dbgetv(dbors,'orid','lat','lon','depth','time','ms');
dbclose(db);
evdir = [num2str(orids,'%03d'),'_',epoch2str(evtimes,'%Y%m%d%H%M'),'/'];

% prep file download info
fprintf('Looking in %s \n',[datadir,evdir])
stafiles = dir([datadir,evdir]); stafiles = stafiles(3:end);
if isempty(stafiles), error('No station mat files for this event'); end
nstas = length(stafiles);
figure(1), clf, set(gcf,'position',[10 10 1400 1000])
hold on

for is = 1:nstas
    load([datadir,evdir,stafiles(is).name])
    ic = find(strcmp(data.chans.component,comp),1,'first');
    if isempty(ic), continue; end
    datZ = detrend(data.dat(:,ic));
    datZ = datZ/abs(max(datZ));
    plot(data.tt-evtimes,datZ/8+data.gcarc)
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
    text(predar(ind,ip),gcs(ind),phs{ip},'FontSize',25,'color','r')
end