% quick and dirty script to make a record section for a given event

orid = 85;

% antelope db details
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascattendb';

% path to top level of directory tree for data
datadir = '/Users/zeilon/Work/CASCADIA/DATA/'; % needs final slash

%% get to work
db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
dbors = dbsubset(dbor,sprintf('orid == %.0f',orid));
[orids,elats,elons,edeps,evtimes,mag] = dbgetv(dbors,'orid','lat','lon','depth','time','ms');
dbclose(db);
evdir = [num2str(orid,'%03d'),'_',epoch2str(evtime,'%Y%m%d%H%M'),'/'];

% prep file download info
stafiles = dir([datadir,evdir]); stafiles = stafiles(3:end);
nstas = length(stafiles);
figure(1), clf, set(gcf,'position',[10 10 1400 1000])
hold on

for is = 1:nstas
    load([datadir,evdir,stafiles(is).name])
    datZ = detrend(data.dat(:,end));
    datZ = datZ/abs(max(datZ));
    plot(data.tt-evtimes,datZ/8+data.gcarc)
%   
end

lms = axis;
gcs = lms(3):0.2:lms(4);
p = zeros(length(gcs),1); s=p;
for ig = 1:length(gcs)
TT = tauptime('depth',edeps,'phases','P,S,PKS,SKS','deg',gcs(ig));
p(ig) = TT(strcmp({TT.phase},'P')).time;
s(ig) = TT(strcmp({TT.phase},'S')).time;
end

plot(p,gcs,'-r',s,gcs,'--r')