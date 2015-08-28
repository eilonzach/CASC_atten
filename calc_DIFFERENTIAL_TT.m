% cycle through events and calculate spectral ratios for all stations
clear all
addpath('matguts')

%% parameters
phase = 'P';
resamprate = 40 ; % new, common sample rate
component = 'Z'; %'Z', 'R', or 'T'
filtfs = 1./[10 0.5]; % [flo fhi] = 1./[Tmax Tmin] in sec
taperx = 0.2;
xcorwind = [-50 100];
xcorlagmax = 3;
acormin = 0.6;


%% directories 
% DIFF TT RESULTS DIRECTORY
difftdir = '~/Documents/MATLAB/CASC_atten/results_DIFF_TTS/'; % needs final slash
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascattendb';
% DATA DIRECTORY (top level)
datadir = '~/Work/CASCADIA/DATA/'; % needs final slash

%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

% GET EVENTS DATA
db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,elats,elons,edeps,evtimes,mags] = dbgetv(dbor,'orid','lat','lon','depth','time','ms');
norids = dbnrecs(dbor);
dbclose(db);


for ie = 85:85 % 44:norids % loop on orids
    fprintf('\n Orid %.0f %s \n\n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
    evdir = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
    if ~exist([datadir,evdir,'_datinfo.mat'],'file'), fprintf('No station mat files for this event\n');continue, end
    load([datadir,evdir,'_datinfo.mat'])
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end
    
    nstas = length(datinfo);

    ofile = [difftdir,'DT_',num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'_',phase];
    
    cp = struct('samprate',resamprate,'prex',-xcorwind(1),'postx',xcorwind(2),...
            'taperx',taperx,'fhi',filtfs(2),'flo',filtfs(1),'npoles',2,'norm',1);

    
    % RESULTS STRUCTURE
    eq = struct('sta',{datinfo.sta}','clean_filt_parms',[],'dT',0,'acor',0,'abs_arrT',0,'pred_arrT',0);
    eq(1).clean_filt_parms = cp;

    wlen = diff(xcorwind)*resamprate*(1+4*taperx); % window length, accounting for taper+padding
    all_datwf  = zeros(wlen,nstas); % each column will be filtered, windowed, tapered data
%     all_tt     = zeros(wlen,nstas); % each column will be times w.r.t. predicted phase arrival from tauP
    
    % LOOP ON STAS
    for is = 1:nstas
        fprintf('Station %s... ',datinfo(is).sta)        % APPLY DATA QUALITY CONDITIONS
        if datinfo(is).crap == 1; 
            fprintf('crap, skipping\n'),continue 
        end
        if ~datinfo(is).rmresp
            fprintf('must have response removed\n'), continue
        end
        if (strcmp(component,'R') || strcmp(component,'T')) && ~datinfo(is).NEZ
            fprintf('must be rotated\n'), continue
        end
        if strcmp(component,'Z') && ~any(strcmp(datinfo(is).chans,'Z'))
            fprintf('has no vertical\n'), continue
        end
        
        % GET STATION'S DATA
        load([datadir,evdir,datinfo(is).sta,'.mat']); % load sta data for this evt
        eq(is).pred_arrT = data.phases(strcmp({data.phases.phase},phase)).artime;
        
        % ROTATE TO ZRT AND GRAB DESIRED COMPONENT
        datZNE = [data.dat(:,strcmp(data.chans.component,'Z')),...
                  data.dat(:,strcmp(data.chans.component,'N')),...
                  data.dat(:,strcmp(data.chans.component,'E'))];
        datZRT = zne2zrt(datZNE,data.seaz);
        
        switch component
            case 'Z', dat0 = datZRT(:,1);
            case 'R', dat0 = datZRT(:,2);
            case 'T', dat0 = datZRT(:,3);
        end
        
        % RESAMPLE
        ttnew = [data.tt(1):1./resamprate:data.tt(end)]';
        dat0 = interp1(data.tt,dat0,ttnew);
        
        % CLEAN DATA
        cp.pretime = data.phases(strcmp({data.phases.phase},phase)).artime-data.tt(1);
        [ datwf,datf,datc,~,ttsf,tts ] = data_clean( dat0,cp );
        
%         figure(1); clf; hold on
%         plot(ttnew-data.phases(strcmp({data.phases.phase},'P')).artime,dat0./std(dat0),'Linewidth',1.5) % raw data
%         plot(tts,datf./std(datf),'r') % filtered + normalised
%         plot(ttsf,datwf,'g') % filtered, windowed, tapered, padded, normalised
%         xlim(xcorwind)

        % PUT INTO MATRIX
        all_datwf(:,is) = datwf(1:wlen);
        fprintf('got data\n')

                
    end % loop on stas
    
    % DELETE ZERO DATA TRACES
    kill = (mean(abs(all_datwf))==0) | isnan(mean(abs(all_datwf)));
    eq(kill) = [];
    all_datwf(:,kill) = [];
    
    % ITERATIVELY XCOR and remove chans with low acor
    acor = zeros(size(eq));
    good = 1:size(eq);
    while any(acor < acormin)
        % DO CROSS CORRELATION
        [dcor,dcstd,dvcstd,acor]=xcortimes(all_datwf(:,good),1./resamprate,xcorwind(1),xcorlagmax,0);
        % DELETE LOW ACOR CHANS
        good(acor<acormin) = [];
    end
    eq = eq(good);
    
    % STORE RESULTS
    for is = 1:length(eq)
        eq(is).dT       = dcor(is);
        eq(is).acor     = acor(is);
        eq(is).abs_arrT = eq(is).pred_arrT + dcor(is);
    end
    
    fprintf('SAVING...')
    % SAVE
    save(ofile,'eq')
    % PUT BACK INTO DATINFO AND DATA structures
    [~,icor,~] = intersect({datinfo.sta},{eq.sta});
    for is = 1:length(icor)
        fprintf('.')
        load([datadir,evdir,datinfo(icor(is)).sta,'.mat'])
        data.phases(strcmp({data.phases.phase},phase)).xtime = ...
              data.phases(strcmp({data.phases.phase},phase)).time + eq(is).dT;
        data.phases(strcmp({data.phases.phase},phase)).xartime = eq(is).abs_arrT;
        data.phases(strcmp({data.phases.phase},phase)).xacor = eq(is).acor;
        save([datadir,evdir,datinfo(icor(is)).sta],'data')
    end
	[datinfo(icor).xcor] = deal(true);
    save([datadir,evdir,'_datinfo'],'datinfo')
    fprintf('\n')
    return
    
    % RECORD SECTION TO CHECK
    
    
    
    
    
    
end