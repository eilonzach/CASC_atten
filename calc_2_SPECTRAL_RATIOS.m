% cycle through events and calculate spectral ratios for all stations
% clear all
close all
cd /Users/zeilon/Documents/MATLAB/CASC_atten
addpath('matguts')

%% parameters
phase = 'S';
component = 'T'; %'Z', 'R', or 'T'
resamprate = 5 ; % new, common sample rate
filtfs = 1./[40 .5]; % [flo fhi] = 1./[Tmax Tmin] in sec
taperx = 0.2;
datwind = [-160 165]; % window of data in eqar structure
specwind = [-5 30];
snrmin = 5;
mavwind = 1; % length of moving average to smooth spectrum (1 for no smoothing)
lofrq = 0.05; % uppermost freq to fit (Hz)
hifrq = 0.25; % uppermost freq to fit (Hz)

overwrite = true;
ifOBSonly = false;
ifusecorZ = true;
ifplot    = false;
ifsave    = true;

%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% DATA DIRECTORY (top level)
datadir = '/Volumes/DATA_mini2/CASCADIA/DATA/'; % needs final slash

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

for ie = 1:270 % 44:norids % loop on orids
%     if  mags(ie)<6.9, continue, end
    tic
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
    if ~isfield(datinfo,'xcor') 
        fprintf('Need to xcor arrival time for this phase and event\n'), continue
    end
	if ~any([datinfo.xcor]==true)
        fprintf('Need to xcor arrival time for this phase and event\n'), continue
	end
    if isfield(datinfo,'dtstar')
        if any([datinfo.dtstar]==true)
            if ~overwrite
                yn = input('delta-tstar already done - overwrite? [y/n] ','s'); 
                if ~strcmp(yn,'y'), fprintf('skipping...\n'), continue, end
            end
        end
    end
    
    %% ------------------ GRAB DATA IN RIGHT FORMAT ------------------
    % prep data structures
    nstas = length(datinfo);
    all_dat0  = zeros(resamprate*diff(datwind),nstas);

    % LOOP ON STAS
    for is = 1:nstas
        fprintf('Station %s... ',datinfo(is).sta)        % APPLY DATA QUALITY CONDITIONS
        if isempty(eqar(is).(['dat' component])),fprintf('no data for this component\n'),continue; end   
        if isempty(eqar(is).dT),fprintf('no xcor for this component\n'),continue; end
        if isnan(eqar(is).dT),fprintf('no xcor for this component\n'),continue; end
        if strcmp(eqar(is).sta,'M08C'), eqar(is).slon = -124.895400; end

        % SHIFT TRACES USING XCOR ARRIVAL TIME 
        % shift so arrival is at time=0
        att = eqar(is).tt-eqar(is).abs_arrT; % shift to since absolute arrival
        ja = (att >= datwind(1)) & (att < datwind(2)); % excerpt times according to datwind
        
        % GRAB DESIRED COMPONENT
        switch component
            case 'Z', all_dat0(:,is) = eqar(is).datZ(ja);
            case 'R', all_dat0(:,is) = eqar(is).datR(ja);
            case 'T', all_dat0(:,is) = eqar(is).datT(ja);
        end
        % grab corrected Z if available and the option is true
        if strcmp(component,'Z') && ifusecorZ && ~isempty(eqar(is).corZ)
            all_dat0(:,is) = eqar(is).corZ(ja);
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
    
    % CALCULATE SNR 
    % from ratio of variance of noise to signal
    snrwf = var(all_datwf)./var(all_datwf_n);
    eqar(1).snr_wf = [];
    eqar = dealto(eqar,'snr_wf',snrwf);

    % ONLY USE GOOD TRACES
    indgd = 1:size(eqar);
    indgd(mean(abs(all_datwf(:,indgd)))==0)     = []; % kill zero traces
    indgd(isnan(mean(abs(all_datwf(:,indgd))))) = []; % kill nan traces
    indgd(snrwf(indgd)<snrmin)                  = []; % kill low snr traces
    if length(indgd) < 2, fprintf('NO GOOD TRACES/ARRIVALS, skip...\n'), continue, end
    
    
    %% ------------------ CALCULATE SPECTRA ------------------
    % WORK OUT SOME KEY VALS
 	wlen = diff(specwind)*resamprate*(1+4*taperx); % window length, accounting for taper+padding
    dt = 1./resamprate;
    ntap=2;
    nft=2^nextpow2(wlen);
    nyq=0.5*resamprate;

    hw = waitbar(0,'PROGRESS THROUGH CALCULATION OF SPECTRA');
    for ig = 1:length(indgd)
        is = indgd(ig);
        
        [specn,~]=pmtm(all_datwf_n(:,is),ntap,nft,resamprate);
        [specs,frq]=pmtm(all_datwf(:,is),ntap,nft,resamprate);
        frq=frq(2:length(frq));
        specn=specn(2:length(specn));
        specs=specs(2:length(specs));

        eqar(is).frq = frq;
% %         %convert power to amplitude, and integrate to displacement
% %         eqar(is).specn=(specn.^0.5)./(2.*pi.*frq);
% %         eqar(is).specs=(specs.^0.5)./(2.*pi.*frq);
        %convert power to amplitude ALREADY IN DISPLACEMENT
        eqar(is).specn=(specn.^0.5);
        eqar(is).specs=(specs.^0.5);
        eqar(is).specss = moving_average(eqar(is).specs,mavwind);

        eqar(is).fcross = frq(max([find(eqar(is).specs<eqar(is).specn,1,'first'),2])-1);
        if isempty(eqar(is).fcross), eqar(is).fcross = 0; end
        waitbar(ig/length(indgd),hw)
    end
    delete(hw)
%     % ONLY USE GOOD TRACES
%     indgd([eqar(indgd).fcross]<hifrq)  	= []; % kill low f-noise crossing traces
    if length(indgd) < 3, fprintf('NOT ENOUGH GOOD TRACES/ARRIVALS, skip...\n'), continue, end
    
 	%% ------------------ CALCULATE DT-STAR FROM REFSTA ------------------
%     refsta = 'WISH';
%     iref = find(strcmp({eqar.sta},refsta));
%     refspecss = eqar(iref).specss;
%     
%     lnR_all = zeros(length(eqar(is).specss),length(indgd));
% 	ind = frq<hifrq;
%     figure(78), clf, set(gcf,'position',[440 0 1000 1400]), hold on
%     for ig = 1:length(indgd)
%         is = indgd(ig);
%         ispecss = eqar(is).specss;
%         lnR = log(ispecss./refspecss);
%         fo = fit(frq(ind),lnR(ind),'poly1');
%         
%         subplot(2,1,1), hold on
%         hr = plot(frq,lnR,'Linewidth',1.5);
%         hrf = plot(fo);
%         xlim([0 hifrq])
%         
%         subplot(2,1,2), hold on
%         hs = plot(ttws,all_datwf(:,is)./max(max(abs(all_datwf(:,indgd))))+eqar(is).gcarc,'Linewidth',1.5);
%         
%         set(hs,'color',colour_get(eqar(is).slon,max([eqar.slon]),min([eqar.slon])))
%         set(hr,'color',colour_get(eqar(is).slon,max([eqar.slon]),min([eqar.slon])))
%         set(hrf,'color',colour_get(eqar(is).slon,max([eqar.slon]),min([eqar.slon])))
%         
%         lnR_all(:,ig) = lnR;
%     end
    
    %% ------------------ CALCULATE DIFFERENTIAL T-STAR ------------------
    specss = zeros(length(eqar(is).specss),length(indgd));
    for ig = 1:length(indgd)
        is = indgd(ig);
        specss(:,ig) = eqar(is).specss;
    end
    
    % CALC DELTA-TSTAR.
    fprintf('Calculate least-squares differential t-star\n')
    [ delta_tstar,cov_dtstar,std_dtstar ] = xspecratio( specss,frq,hifrq,lofrq,1,ifplot );
    
    
    %% ---------------------- STORE RESULTS -----------------------
    % STORE RESULTS
    fprintf('Recording results in arrival structure...')
    % prep eqar to receive new fields
    eqar(1).dtstar = []; eqar(1).std_dtstar = []; eqar(1).par_dtstar_specR = [];
    par_dtstar = struct('comp',component,'filtfs',filtfs,'window',specwind,...
                          'taperx',taperx,'mavwind',mavwind,'hifrq',hifrq,'lofrq',lofrq,'snrmin',snrmin);
                      
    eqar(indgd) =  dealto(eqar(indgd),'dtstar',delta_tstar);
    eqar(indgd) =  dealto(eqar(indgd),'std_dtstar',std_dtstar);
    eqar(indgd) =  dealto(eqar(indgd),'par_dtstar_specR',par_dtstar);

	indbd = setdiff(1:length(eqar),indgd);
    if ~isempty(indbd)
    eqar(indbd) =  dealto(eqar(indbd),'dtstar',nan);
    end
    
    %% -------------------------- PLOTS ---------------------------
    if ifplot
    %% look at maximum frequencies above noise
    figure(87), clf, hold on
    fcross = [eqar(indgd).fcross]';
    hist(fcross,100); xlim([0 2])
    line(hifrq*[1 1],[0 max(hist(fcross))],'Color','r','LineStyle','--','LineWidth',2)
    title('maximum freq where signal is above noise')

    %% spectral plot
    figure(1), clf, set(gcf,'position',[440 0 1000 1400]), hold on
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
        hps = semilogx(frq(1:n4),log10(eqar(is).specss(1:n4)),'b','LineWidth',1.5); % << here using smoothed spectrum 
        hpn = semilogx(frq(1:n4),log10(eqar(is).specn(1:n4)),'k','LineWidth',1.5);
        set(hps,'color',colour_get(eqar(is).slon,max([eqar.slon]),min([eqar.slon])))
        % xlim(ax(1),[0 0.8*nyq]);
        xlabel('Hz'); 
        ylabel('amplitude, nm/Hz')
        % a=axis(ax(1));
        % % ylim_max=10^(floor(log10(max(specs(1:n4))))+1);
        % % ylim(ax(1),ylim_max.*[1.e-4 1])
        xlim([0 hifrq]);
        ylim([-9,-2]);
    end % loop on good stas

    plot_ATTEN_TandF_domain( eqar(indgd) )
    end % ifplot
    
	%% -------------------------- SAVE ---------------------------
    if ifsave
    % SAVE
    save(arfile,'eqar')
    fprintf(' saved\n')
    % RECORD DTSTARCALC IN DATINFO 
    [datinfo.dtstar] = deal(false);
	[datinfo(indgd).dtstar] = deal(true);
    save(datinfofile,'datinfo')
    
    fprintf(' STA  CHAN  NEZ  resp  tilt  comp  xcor  tstar\n')
	for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f     %1.0f     %1.0f     %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).rmtilt,datinfo(is).rmcomp,datinfo(is).xcor,datinfo(is).dtstar); end
    end

    
    toc  
end % loop on orids

% results_PARSE