% cycle through events and calculate differential travel times for all stations
clear all
close all
addpath('matguts')

%% parameters
phase = 'P';
resamprate = 40 ; % new, common sample rate
component = 'Z'; %'Z', 'R', or 'T'
filtfs = 1./[4 0.1]; % [flo fhi] = 1./[Tmax Tmin] in sec
taperx = 0.2;
datwind = [-200 200]; % window of data in eqar structure
xcorwind = [-40 20];
xcorlagmax = 3;
acormin = 0.6;


%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascattendb';
% DATA DIRECTORY (top level)
datadir = '~/Work/CASCADIA/DATA/'; % needs final slash

%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

cp = struct('samprate',resamprate,'pretime',-datwind(1),'prex',-xcorwind(1),'postx',xcorwind(2),...
            'taperx',taperx,'fhi',filtfs(2),'flo',filtfs(1),'npoles',2,'norm',1);

% GET EVENTS DATA
db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,elats,elons,edeps,evtimes,mags] = dbgetv(dbor,'orid','lat','lon','depth','time','ms');
norids = dbnrecs(dbor);
dbclose(db);


for ie = 85:85 % 44:norids % loop on orids
    tic
    fprintf('\n Orid %.0f %s \n\n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
    % name files and directories
    evdir    = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
    datinfofile = [datadir,evdir,'_datinfo'];
    arfile   = [datadir,evdir,'_EQAR_',phase];
   
    % check files exist
    if ~exist([datinfofile,'.mat'],'file'), fprintf('No station mat files for this event\n');continue, end
    if ~exist([arfile,'.mat'],'file'), fprintf('No arfile for this event + phase\n');continue, end
    
    % load files
    load(datinfofile) % loads datinfo stucture
    load(arfile)      % loads eqar structure
    
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end
    nstas = length(datinfo);
    
    %% ------------------ GRAB DATA IN RIGHT FORMAT ------------------
    % prep data structures
    wlen = diff(xcorwind)*resamprate*(1+4*taperx); % window length, accounting for taper+padding
    all_dat0  = zeros(resamprate*diff(datwind),nstas);
    
    % LOOP ON STAS
    for is = 1:nstas
        fprintf('Station %s... ',datinfo(is).sta)        % APPLY DATA QUALITY CONDITIONS
        if isempty(eqar(is).(['dat' component])),fprintf('no data for this component\n'),continue; end   
        
        % GRAB DESIRED COMPONENT
        switch component
            case 'Z', all_dat0(:,is) = eqar(is).datZ(:);
            case 'R', all_dat0(:,is) = eqar(is).datR(:);
            case 'T', all_dat0(:,is) = eqar(is).datT(:);
        end
        fprintf('got data\n')
    end % loop on stas
    % CLEAN DATA
    [ all_datwf,all_datf,all_datwc,all_datc,~,ttws,tts ] = data_clean( all_dat0,cp );

    % MARK ZERO DATA TRACES
    indgd = 1:size(eqar);
    indgd(mean(abs(all_datwf(:,indgd)))==0)     = []; % kill zero traces
    indgd(isnan(mean(abs(all_datwf(:,indgd))))) = []; % kill nan traces
    if length(indgd) < 2, fprintf('NO GOOD TRACES/ARRIVALS, skip...\n'), continue, end
    
    %% ------------------ ITERATIVELY XCOR ------------------
    % removing chans with low acor
    fprintf('Cross-correlating + iteratively removing low acor traces\n')
    acor = zeros(size(eqar));
    
    while any(acor < acormin)
        % DO CROSS CORRELATION
        [dcor,dcstd,dvcstd,acor]=xcortimes(all_datwf(:,indgd),1./resamprate,xcorwind(1),xcorlagmax,1);
        % DELETE LOW ACOR CHANS
        indgd(acor<acormin) = []; % kill low acor traces
        if length(indgd) < 2, break, end
    end
	if length(indgd) < 2, fprintf('NO GOOD TRACES/ARRIVALS, skip...\n'), continue, end
    
    %% ------------------ PLOTTING ------------------
    % RECORD SECTION TO CHECK
    figure(54), clf, set(gcf,'position',[10 600 1200 400]), hold on
    normf = max(max(abs(all_datwf(:,indgd))));
    gcarcs = [eqar.gcarc];
    for ig = 1:length(indgd)
        is = indgd(ig);
        plot(ttws                   ,all_datwf(:,is)/normf/2 + gcarcs(is),'k')
        plot(ttws-dcor(ig),all_datwf(:,is)/normf/2 + gcarcs(is),'b','LineWidth',1.5)
        text(ttws(end)+0.5,gcarcs(is),datinfo(is).sta)
    end
	% STACK raw traces to get overall arrtime
    stk = zeros(wlen+1,1);
    for ig = 1:length(indgd)
        is = indgd(ig);
        nshift = round(dcor(ig)*resamprate);
        % can only do this trick because padded with zeros
        datstk = [zeros(-2*nshift,1);all_datwf(abs(nshift)+1:end-abs(nshift),is);zeros(2*nshift,1)];
        stk = stk + datstk;
    end
    figure(54), hold on
    axlim = axis;
    plot(ttws(1:wlen+1),stk./max(abs(stk)) + axlim(3),'r','LineWidth',1.5)
	text(ttws(end)+0.5,axlim(3),'*STACK*')

    % PICK phase arrival
    axlim = axis;
    conf = 0; hl = [];
    while conf==0
        title('PICK PHASE ARRIVAL TIME (CLICK OUTSIDE PLOT TO ABORT)','FontSize',18)
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
    
%% ---------------------- STORE RESULTS -----------------------
    % STORE RESULTS
    fprintf('Recording results in arrival structure...')
    % prep eqar to receive new fields
    eqar(1).dT = []; eqar(1).acor = []; eqar(1).abs_arrT = []; eqar(1).par_dT = [];

    eqar(indgd) =  dealto(eqar(indgd),'dT',dcor);
    eqar(indgd) =  dealto(eqar(indgd),'acor',acor);
    eqar(indgd) =  dealto(eqar(indgd),'abs_arrT',[eqar(indgd).pred_arrT]' + dcor + t0(1));
    eqar(indgd) =  dealto(eqar(indgd),'par_dT',cp);

    pause(.1)

    
      
    % SAVE
    save(arfile,'eqar')
    fprintf(' saved\n')
    % RECORD XCOR IN DATINFO 
	[datinfo(indgd).xcor] = deal(true);
    save([datadir,evdir,'_datinfo'],'datinfo')
    
    toc  
    return
end