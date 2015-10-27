% clear all
close all
cd('/Users/zeilon/Documents/MATLAB/CASC_atten')
addpath('/Users/zeilon/Documents/MATLAB/seis_tools') %make sure this dir's functions get used
addpath('matguts') %make sure this dir's functions get used

respdir = '/Users/zeilon/Work/CASCADIA/CAdb/response_BarclayPZ/';
dbpath = '/Users/zeilon/Work/CASCADIA/CAdb/cascattendb';
dbrespdir = '/Users/zeilon/Work/CASCADIA/CAdb/response_new/';
seedrespdir = '~/Work/CASCADIA/CAdb/RESP_NW_STA_CHA/';

fmin = 1e-4;
fmax = 1e3;

%% e.g. station
% LDEO trawl: FN07A (ins 1,2)
% Scripps Abalone: J35A (ins 3,4)
% WHOI: FS05B (ins 5,6), G04B (ind 7,8)
sta = 'FN07A'; 
chan = 'HHZ';
net = '7D';


if strcmp(chan(2:3),'DH'),extra_zp=0;typstr='';    end
if strcmp(chan(2),'H'),   extra_zp=1;typstr='VEL'; end
% extra_zp = 0;

% IRIS respfile
db = dbopen(dbpath,'r');
dbsn = dblookup_table(db,'sensor');
dbin = dblookup_table(db,'instrument');
dbsns = dbsubset(dbsn,sprintf('sta=="%s" && chan=="%s"',sta,chan));
if dbnrecs(dbsns)==0, fprintf('That sta+chan doesn''t exist\n'); dbclose(db); return; end
dbj1 = dbjoin(dbsns,dbin);
[rfile] = dbgetv(dbj1,'dfile');
dbclose(db);
respfile_IRIS = [dbrespdir,rfile];

% Barclay respfile
rfile = dir([respdir,'*',sta,'_',chan,'*']);
respfile_AB = [respdir,rfile.name];

% FULL respfile (from dataless seed)
respfilestr = ['RESP.*.',sta,'..',chan ];
respfile = dir([seedrespdir,respfilestr]);

%% plot IRIS response
[zz,pp,gain] = read_sac_pole_zero_haj(respfile_IRIS);
plot_resp_pz( pp,zz,gain, fmin, fmax,extra_zp )
subplot(211),title(sprintf('IRIS %s response for Station %s, Channel %s',typstr,sta,chan),'FontSize',18)
clone_figure(33,34)

%% plot LDEO response
[zz,pp,gain] = read_sac_pole_zero_haj(respfile_AB);
plot_resp_pz(pp,zz,gain, fmin, fmax,extra_zp )
subplot(211),title(sprintf('LDEO %s response for Station %s, Channel %s',typstr,sta,chan),'FontSize',18)


%% plot FULL response
[ resp,faxis ] = mkresp_evalresp( seedrespdir,net,sta,chan,1e4,1e7,'2011 300' );
plot_resp( resp,faxis,extra_zp )
subplot(211),title(sprintf('FULL %s response for Station %s, Channel %s',typstr,sta,chan),'FontSize',18)

return
%% e.g. event
orid = 150;

db = dbopen(dbpath,'r');
dbor = dblookup_table(db,'origin');
dbor.record = dbfind(dbor,sprintf('orid == %.0f',orid));
evtime = dbgetv(dbor,'time');
dbclose(db);
evdir = ['/Volumes/DATA/CASCADIA/DATA/',num2str(orid),'_',epoch2str(evtime,'%Y%m%d%H%M'),'/'];
if ~exist([evdir,sta,'.mat'],'file'), error('No data for this sta+orid'); end

load([evdir,sta,'.mat']);
ic = strcmp(data.raw.chans.name,chan);

[zz,pp,gain] = read_sac_pole_zero_haj(respfile_IRIS);
odat = rm_resp(data.raw.dat(:,ic),zz,pp,gain,data.samprate,extra_zp,1);
subplot(3,2,5:6),title(sprintf('IRIS %s response for orid %.0f, Station %s, Channel %s',typstr,orid,sta,chan),'FontSize',18)
clone_figure(39,40)


[zz,pp,gain] = read_sac_pole_zero_haj(respfile_AB);
odat = rm_resp(data.raw.dat(:,ic),zz,pp,gain,data.samprate,extra_zp,1);
subplot(3,2,5:6),title(sprintf('LDEO %s response for orid %.0f, Station %s, Channel %s',typstr,orid,sta,chan),'FontSize',18)

