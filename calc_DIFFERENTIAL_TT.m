% cycle through events and calculate differential travel times for all stations
clear all
close all
addpath('matguts')

%% parameters
phase = 'P';
component = 'Z'; %'Z', 'R', or 'T'
resamprate = 40 ; % new, common sample rate
% Do filtfs [12 1] for S, and [5 1] for P
filtfs = 1./[5 1]; % [flo fhi] = 1./[Tmax Tmin] in sec 
taperx = 0.2;
datwind = [-200 200]; % window of data in eqar structure
prelimwind = [-40 40];
% prelimwind = [-25 15];
xcorlagmax = 6;
acormin = 0.65;
overwrite = true;

ifOBSonly = true;

manualkillstas = {''};

%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% DATA DIRECTORY (top level)
datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash

%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%


% GET EVENTS DATA
db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,elats,elons,edeps,evtimes,mags] = dbgetv(dbor,'orid','lat','lon','depth','time','ms');
norids = dbnrecs(dbor);
dbclose(db);

obsstr = ''; if ifOBSonly, obsstr = 'OBS_'; end

for ie = 269:norids% 44:norids % loop on orids % got to 475 on SKS,R, % got to 190 on S,T
%     if  mags(ie)<6.9, continue, end
    tic
    close all
    fprintf('\n Orid %.0f %s \n\n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
    % name files and directories
    evdir       = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
    datinfofile = [datadir,evdir,'_datinfo_',obsstr,phase];
    arfile      = [datadir,evdir,'_EQAR_',obsstr,phase,'_',component];
   
    % check files exist
    if ~exist([datinfofile,'.mat'],'file'), fprintf('No station mat files for this event\n');continue, end
    if ~exist([arfile,'.mat'],'file'), fprintf('No arfile for this event + phase\n');continue, end
    
    % load files
    load(datinfofile) % loads datinfo stucture
    load(arfile)      % loads eqar structure
    
    % options to skip
    if isempty(datinfo), fprintf('No station mat files for this event\n'); continue, end
    if isfield(datinfo,'xcor')
        if any([datinfo.xcor]==true)
            if ~overwrite
                yn = input('Xcor already done - overwrite? [y/n] ','s'); 
                if ~strcmp(yn,'y'), fprintf('skipping...\n'), continue, end
            end
        end
    end
    if all([eqar.pred_arrT]==0), fprintf('No %s arrivals for this event\n',phase), continue, end
    
    if strcmp(phase,'SKS') && mean([eqar.gcarc])<85, continue, end
	if strcmp(phase,'SKP') && mean([eqar.gcarc])<129, continue, end

    %% ------------------ GRAB DATA IN RIGHT FORMAT ------------------
    % prep data structures
    nstas = length(datinfo);
    wlen = diff(prelimwind)*resamprate*(1+4*taperx); % window length, accounting for taper+padding
    all_dat0  = zeros(resamprate*diff(datwind),nstas);
    
    % LOOP ON STAS
    for is = 1:nstas
        fprintf('Station %s... ',datinfo(is).sta)        % APPLY DATA QUALITY CONDITIONS
        if isempty(eqar(is).(['dat' component])),fprintf('no data for this component\n'),continue; end   
        if strcmp(eqar(is).sta,'M08C'), eqar(is).slon = -124.895400; end
        
        % GRAB DESIRED COMPONENT
        switch component
            case 'Z', all_dat0(:,is) = eqar(is).datZ(:);
            case 'R', all_dat0(:,is) = eqar(is).datR(:);
            case 'T', all_dat0(:,is) = eqar(is).datT(:);
        end
        fprintf('got data\n')
    end % loop on stas
    
    % CLEAN DATA
    cp = struct('samprate',resamprate,'pretime',-datwind(1),...
            'prex',-prelimwind(1),'postx',prelimwind(2),...
            'taperx',taperx,'fhi',filtfs(2),'flo',filtfs(1),'npoles',2,'norm',1);
    [ all_datwf,all_datf,all_datwc,all_datc,~,ttws,tts ] = data_clean( all_dat0,cp );

    % MARK ZERO DATA TRACES
    indgd = 1:size(eqar);
    indgd(mean(abs(all_dat0(:,indgd)))<1e-16)     = []; % kill zero traces
    indgd(mean(abs(all_datwf(:,indgd)))==0)     = []; % kill zero traces
    indgd(isnan(mean(abs(all_datwf(:,indgd))))) = []; % kill nan traces
    for ik = 1:length(manualkillstas),indgd(strcmp({eqar(indgd).sta},manualkillstas{ik})) = []; end % kill bad stas traces
    
    if strcmp(phase,'SKS') ,indgd([eqar(indgd).gcarc]<86) = []; end % kill close traces
    
    if length(indgd) < 2, fprintf('NO GOOD TRACES/ARRIVALS, skip...\n'), continue, end
    
    fprintf('Orid %.0f %s \n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
	%% tauptime to help picking
    tauptime('d',mean([eqar.gcarc]),'z',edeps(ie))
  
%% ------------------ PICK XCOR WINDOW, RECLEAN ------------------
    figure(53), clf, set(gcf,'position',[700 100 1200 800]), hold on
    normf = max(max(abs(all_datwf(:,indgd))));
    for ig = 1:length(indgd)
        is = indgd(ig);
        plot(ttws,all_datwf(:,is)/normf/2 + eqar(is).gcarc,'--k','LineWidth',.5)
        text(prelimwind(2)+1,eqar(is).gcarc,datinfo(is).sta)
    end
    axis([prelimwind min([eqar(indgd).gcarc])-0.5 max([eqar(indgd).gcarc])+1 ])
    title('PICK XCOR WINDOW','FontSize',18)
    xcorwind = ginput(2); 
    xcorwind = round_level(xcorwind(:,1),1./resamprate);
    
    % CLEAN DATA
    cp = struct('samprate',resamprate,'pretime',-datwind(1),...
            'prex',-xcorwind(1),'postx',xcorwind(2),...
            'taperx',taperx,'fhi',filtfs(2),'flo',filtfs(1),'npoles',2,'norm',1);
    [ xcor_datwf,~,~,~,~,xcor_ttws,~ ] = data_clean( all_dat0,cp );

    % MARK ZERO DATA TRACES
    indgd(mean(abs(all_dat0(:,indgd)))<1e-16)     = []; % kill zero traces
    indgd(mean(abs(xcor_datwf(:,indgd)))==0)     = []; % kill zero traces
    indgd(isnan(mean(abs(xcor_datwf(:,indgd))))) = []; % kill nan traces
    if length(indgd) < 2
        fprintf('NO GOOD TRACES/ARRIVALS, skip...\n')
        [datinfo.xcor] = deal(false);
        save(datinfofile,'datinfo');
        continue
    end
    
    %% ------------------ ITERATIVELY XCOR ------------------
    % removing chans with low acor
    fprintf('Cross-correlating + iteratively removing low acor traces\n')
    acor = zeros(size(eqar));
    
    acor_all = zeros(size(indgd));
    
    iter = 1;
    while any(acor < acormin)
        % DO CROSS CORRELATION
        [dcor,dcstd,dvcstd,acor]=xcortimes(xcor_datwf(:,indgd),1./resamprate,xcorwind(1),xcorlagmax,1);
        
        acor_all(indgd) = acor;
        
        % DELETE LOW ACOR CHANS
        if iter == 1, 
            indgd(acor < 0.1) = []; % first kill really bad traces
        elseif iter > 1
            indgd(acor<acormin) = []; % then kill low acor traces
        end           
        if length(indgd) < 2, break, end
        iter=iter+1;
    end
	if length(indgd) < 2
        fprintf('NO GOOD TRACES/ARRIVALS, skip...\n')
        [datinfo.xcor] = deal(false);
        save(datinfofile,'datinfo');
        continue
    end
    
    %% ------------------ PLOTTING ------------------
    % RECORD SECTION TO CHECK
    figure(54), clf, set(gcf,'position',[10 100 1200 800]), hold on
    normf = max(max(abs(all_datwf(:,indgd))));
    normfx = max(max(abs(xcor_datwf(:,indgd))));
    gcarcs = [eqar.gcarc];
    
    kk = 0;
    for ig = 1:length(indgd)
        is = indgd(ig);
        kk = kk+1;
%         plot(ttws                   ,all_datwf(:,is)/normf/2 + gcarcs(is),'k')
        plot(ttws-dcor(ig),all_datwf(:,is)/normf + kk,'b','LineWidth',1.5)
        plot(xcor_ttws-dcor(ig),xcor_datwf(:,is)/normfx + kk,'g','LineWidth',1.5)
        text(ttws(end)+1,kk,datinfo(is).sta)
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
	xlabel(sprintf('Orid %.0f,  gcarc ~%.1f,  seaz ~%.1f',ie,mean([eqar.gcarc]),mean([eqar.seaz])),'Fontsize',22)
    ylim([-1 1+length(indgd)])
   
    % PICK phase arrival
    axlim = axis;
    conf = 0; hl = [];
    while conf==0
        title('PICK PHASE ARRIVAL TIME (CLICK OUTSIDE PLOT TO SKIP)','FontSize',18)
        delete(hl);
        t0 = ginput(1);
        if t0(1) > axlim(2) || t0(1) < axlim(1)
            fprintf('Skipping event on user request...\n')
            conf=nan;
        else
            hold on 
            hl = plot([t0(1),t0(1)],[axlim(3),axlim(4)],'--k');
            title('CLICK IN PLOT TO CONFIRM','FontSize',18)
            pause(0.01)
            conf = ginput(1);
            if conf(1) > axlim(2) || conf(1) < axlim(1), conf = 0; end
        end
    end
    if isnan(conf), 
        [datinfo.xcor] = deal(false);
        save(datinfofile,'datinfo');
        continue
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

	indbd = setdiff(1:length(eqar),indgd);
    if ~isempty(indbd)
        eqar(indbd) =  dealto(eqar(indbd),'dT',nan);
    end

    
    pause(.1)

    
      
    % SAVE
    save(arfile,'eqar')
    fprintf(' saved\n')
    % RECORD XCOR IN DATINFO 
    [datinfo.xcor] = deal(false);
    [datinfo(indgd).xcor] = deal(true);
    save(datinfofile,'datinfo')
    
    fprintf(' STA  CHAN  NEZ  resp  tilt  comp  xcor\n')
	for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f     %1.0f     %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).rmtilt,datinfo(is).rmcomp,datinfo(is).xcor); end
    toc  

end