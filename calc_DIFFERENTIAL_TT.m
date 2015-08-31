% cycle through events and calculate spectral ratios for all stations
clear all
% close all
addpath('matguts')

%% parameters
phase = 'P';
resamprate = 40 ; % new, common sample rate
component = 'Z'; %'Z', 'R', or 'T'
filtfs = 1./[4 0.1]; % [flo fhi] = 1./[Tmax Tmin] in sec
taperx = 0.2;
xcorwind = [-40 20];
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
    tic
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
    all_datwc  = zeros(wlen,nstas); % each column will be clean, windowed, tapered data
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
        [ datwf,datf,datwc,datc,~,ttws,tts ] = data_clean( dat0,cp );
        
%         figure(1); clf; hold on
%         plot(ttnew-data.phases(strcmp({data.phases.phase},'P')).artime,dat0./std(dat0),'Linewidth',1.5) % raw data
%         plot(ttws,datwc./std(datwc),'k','Linewidth',1.5) % filtered, windowed, tapered, padded, normalised
%         plot(tts,datf./std(datf),'r') % filtered + normalised
%         plot(ttws,datwf,'g') % filtered, windowed, tapered, padded, normalised
%         xlim(xcorwind)

        % PUT INTO MATRIX
        all_datwf(:,is) = datwf(1:wlen);
        all_datwc(:,is) = datwc(1:wlen);
        fprintf('got data\n')

                
    end % loop on stas
    
    % MARK ZERO DATA TRACES
    indgd = 1:size(eq);
    indgd(mean(abs(all_datwf(:,indgd)))==0)     = []; % kill zero traces
    indgd(isnan(mean(abs(all_datwf(:,indgd))))) = []; % kill nan traces
    
    % ITERATIVELY XCOR and remove chans with low acor
    fprintf('Cross-correlating + iteratively removing low acor traces\n')
    acor = zeros(size(eq));
    
    while any(acor < acormin)
        % DO CROSS CORRELATION
        [dcor,dcstd,dvcstd,acor]=xcortimes(all_datwf(:,indgd),1./resamprate,xcorwind(1),xcorlagmax,1);
        % DELETE LOW ACOR CHANS
        indgd(acor<acormin) = []; % kill low acor traces
    end
    
	% DELETE NOT GOOD TRACES
    eq = eq(indgd);
    % STORE RESULTS
    for is = 1:length(eq)
        eq(is).dT       = dcor(is);
        eq(is).acor     = acor(is);
        eq(is).abs_arrT = eq(is).pred_arrT + dcor(is);
    end
    pause(1)

    
    % RECORD SECTION TO CHECK
    figure(54), clf, set(gcf,'position',[10 600 1200 400]), hold on
    for is = 1:length(indgd)
        load([datadir,evdir,datinfo(indgd(is)).sta,'.mat'])
        plot(ttws(1:wlen) ,all_datwf(:,indgd(is))/30+data.gcarc,'k')
        plot(ttws(1:wlen)-eq(is).dT ,all_datwf(:,indgd(is))/30+data.gcarc,'b','LineWidth',1.5)
        text(ttws(end)+0.5,data.gcarc,data.station.name)
    end

	% STACK raw traces to get overall arrtime
    stk = zeros(wlen,1);
    for is = 1:indgd
        nshift = round(dcor(is)*resamprate);
        % can only do this trick because padded with zeros
        datstk = [zeros(-2*nshift,1);all_datwf(abs(nshift)+1:end-abs(nshift),indgd(is));zeros(2*nshift,1)];
        stk = stk + datstk;
    end
    figure(54), hold on
    axlim = axis;
    plot(ttws(1:wlen),stk./30 + axlim(3),'r','LineWidth',1.5)

    % PICK phase arrival
    axlim = axis;
    conf = 0; hl = [];
    while conf==0
        title('PICK PHASE ARRIVAL TIME (CLICK OUTSIDE PLOT TO ABORT','FontSize',18)
        delete(hl);
        t0 = ginput(1);
        if t0(1) > axlim(2) || t0(1) < axlim(1)
            error('Aborting xcorr on user request...')
        else
            hold on 
            hl = plot([t0(1),t0(1)],[axlim(3),axlim(4)],'--k');
            title('CLICK IN PLOT TO CONFIRM','FontSize',18)
            pause(0.01)
            conf = ginput(1);
            if conf(1) > axlim(2) || conf(1) < axlim(1), conf = 0; end
        end
    end
    % save into eq
    temp = num2cell([eq.abs_arrT] + t0(1));
    [eq.abs_arrT] = deal(temp{:});    

    hw = waitbar(0,'SAVING...');
    % SAVE
    save(ofile,'eq')
    % PUT BACK INTO DATINFO AND DATA structures
    for is = 1:length(indgd)
        load([datadir,evdir,datinfo(indgd(is)).sta,'.mat'])
        ph = data.phases; ip = strcmp({ph.phase},phase);
        ph(ip).xtime = ph(ip).time + eq(is).dT;
        ph(ip).xartime = eq(is).abs_arrT;
        ph(ip).xacor = eq(is).acor;
        data.phases = ph;
        save([datadir,evdir,data.station.name],'data')
        waitbar(is/length(indgd),hw)
    end
	[datinfo(indgd).xcor] = deal(true);
    save([datadir,evdir,'_datinfo'],'datinfo')
    delete(hw)
    
    toc
    return
    
    
    
    
    
    
end