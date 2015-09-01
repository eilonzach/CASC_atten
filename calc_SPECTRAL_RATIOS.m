% cycle through events and calculate spectral ratios for all stations
clear all
close all
addpath('matguts')

%% parameters
phase = 'P';
resamprate = 40 ; % new, common sample rate
component = 'Z'; %'Z', 'R', or 'T'
filtfs = 1./[40 .1]; % [flo fhi] = 1./[Tmax Tmin] in sec
taperx = 0.2;
datwind = [-180 175]; % window of data in eqar structure
specwind = [-10 30];
snrmin = 3;
mavwind = 1; % length of moving average to smooth spectrum (1 for no smoothing)
hifrq = 0.4; % uppermost freq to fit (Hz)

%% directories 
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
    if ~any([datinfo.xcor]), fprintf('Need to xcor arrival time for this phase and event\n');continue, end
    nstas = length(datinfo);
    
    %% ------------------ GRAB DATA IN RIGHT FORMAT ------------------
    % prep data structures
    all_dat0  = zeros(resamprate*diff(datwind),nstas);
    
    % LOOP ON STAS
    for is = 1:nstas
        fprintf('Station %s... ',datinfo(is).sta)        % APPLY DATA QUALITY CONDITIONS
        if isempty(eqar(is).(['dat' component])),fprintf('no data for this component\n'),continue; end   
        if isempty(eqar(is).dT),fprintf('no xcor for this component\n'),continue; end
        
        att = eqar(is).tt-eqar(is).abs_arrT; % shift to since absolute arrival
        ja = (att >= datwind(1)) & (att < datwind(2)); % excerpt times according to datwind
        
        % GRAB DESIRED COMPONENT
        switch component
            case 'Z', all_dat0(:,is) = eqar(is).datZ(ja);
            case 'R', all_dat0(:,is) = eqar(is).datR(ja);
            case 'T', all_dat0(:,is) = eqar(is).datT(ja);
        end
        fprintf('got data\n')
    end % loop on stas
    % CLEAN DATA for signal and noise
    cp = struct('samprate',resamprate,'pretime',-datwind(1),'prex',-specwind(1),'postx',specwind(2),...
            'taperx',taperx,'fhi',filtfs(2),'flo',filtfs(1),'npoles',2,'norm',0);
    cp_n = struct('samprate',resamprate,'pretime',-datwind(1),'prex',-specwind(1)+diff(specwind),'postx',specwind(2)-diff(specwind),...
            'taperx',taperx,'fhi',filtfs(2),'flo',filtfs(1),'npoles',2,'norm',0);

    
    [ all_datwf,all_datf,all_datwc,all_datc,~,ttws,tts ] = data_clean( all_dat0,cp );
    [ all_datwf_n,all_datf_n,all_datwc_n,all_datc_n,~,ttws_n,tts_n ] = data_clean( all_dat0,cp_n );
    
    
    % MARK ZERO DATA TRACES
    indgd = 1:size(eqar);
    indgd(mean(abs(all_datwf(:,indgd)))==0)     = []; % kill zero traces
    indgd(isnan(mean(abs(all_datwf(:,indgd))))) = []; % kill nan traces
    if length(indgd) < 2, fprintf('NO GOOD TRACES/ARRIVALS, skip...\n'), continue, end
    
	wlen = diff(specwind)*resamprate*(1+4*taperx); % window length, accounting for taper+padding
    
	snr = var(all_datwf)./var(all_datwf_n);
    
    % WORK OUT SOME KEY VALS
    dt = 1./resamprate;
    ntap=2;
    nft=2^nextpow2(wlen);
    nyq=0.5./dt;

    hw = waitbar(0,'PROGRESS THROUGH SPECTRA CALC.');
    for ig = 1:length(indgd)
        is = indgd(ig);
        
        [specn,~]=pmtm(all_datwf_n(:,is),ntap,nft,resamprate);
        [specs,frq]=pmtm(all_datwf(:,is),ntap,nft,resamprate);
        frq=frq(2:length(frq));
        specn=specn(2:length(specn));
        specs=specs(2:length(specs));

        eqar(is).frq = frq;
        %convert power to amplitude, and integrate to displacement
        eqar(is).specn=(specn.^0.5)./(2.*pi.*frq);
        eqar(is).specs=(specs.^0.5)./(2.*pi.*frq);

        eqar(is).specss = moving_average(eqar(is).specs,mavwind);

        eqar(is).fmax = frq(max([find(eqar(is).specs<eqar(is).specn,1,'first'),2])-1);
        if isempty(eqar(is).fmax), eqar(is).fmax = 0; end
        waitbar(ig/length(indgd),hw)
    end
    delete(hw)
    
    return
    %% look at maximum frequencies above noise
    figure(87), clf, hold on
    fmaxs = [eqar.fmax]';
    hist(fmaxs,100); xlim([0 2])
    line(hifrq*[1 1],[0 max(hist(fmaxs))],'Color','r','LineStyle','--','LineWidth',2)
    title('maximum freq where signal is above noise')

    %% spectral plot
    figure(1), clf
    for ig = 1:length(indgd)
        is = indgd(ig);
	
        % plot data series
        subplot(212), hold on
        hpd = plot(ttws_n,all_datwf_n(:,is),'r',ttws,all_datwf(:,is),'g');
        set(hpd,'color',colour_get(eqar(is).slon,max([eqar.slon]),min([eqar.slon])))
        xlabel('time from pick, s')
        
        % plot spectrum
        subplot(211), hold on
        n4=length(find(frq<=0.8*nyq));
        hps = plot(frq(1:n4),log10(eqar(is).specss(1:n4)),'b','LineWidth',1.5); % << here using smoothed spectrum 
        hpn = plot(frq(1:n4),log10(eqar(is).specn(1:n4)),'k','LineWidth',1.5);
        set(hps,'color',colour_get(eqar(is).slon,max([eqar.slon]),min([eqar.slon])))
        % xlim(ax(1),[0 0.8*nyq]);
        xlabel('Hz'); 
        ylabel('amplitude, nm/Hz')
        % a=axis(ax(1));
        % % ylim_max=10^(floor(log10(max(specs(1:n4))))+1);
        % % ylim(ax(1),ylim_max.*[1.e-4 1])
        xlim([0 hifrq]);
        ylim([-9,-2]);
    end

    
    
    %     % STORE RESULTS
%     fprintf('Recording results in arrival structure...')
%     for is = 1:length(indgd)
%         eqar(indgd(is)).dT       = dcor(is);
%         eqar(indgd(is)).acor     = acor(is);
%         eqar(indgd(is)).abs_arrT = eqar(is).pred_arrT + dcor(is) + t0(1);
%         eqar(indgd(is)).dT_filt_parms = cp;
%     end
%     pause(.1)
%     
%       
%     % SAVE
%     save(arfile,'eqar')
%     fprintf(' saved\n')
%     % RECORD XCOR IN DATINFO 
% 	[datinfo(indgd).xcor] = deal(true);
%     save([datadir,evdir,'_datinfo'],'datinfo')
%     
    toc  
    return
end