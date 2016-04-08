% cycle through events and calculate spectral ratios for all stations
% pause
clear all
close all
cd /Users/zeilon/Documents/MATLAB/CASC_atten
addpath('matguts')

%% parameters
phase = 'S';
component = 'T'; %'Z', 'R', or 'T'
resamprate = 40 ; % new, common sample rate
filtfs = 1./[40 .5]; % [flo fhi] = 1./[Tmax Tmin] in sec
taperx = 0.2;
datwind = [-160 165]; % window of data in eqar structure
specwind = [-5 30];
snrmin = 10;

test_alphas = [0];

overwrite = true;
ifOBSonly = false;
ifusecorZ = false;
ifplot    = true;
ifsave    = true;
ifplotonly = false;

% combspec parms
parms.comb.Tmin = 1;
parms.comb.Tmax = 20;
parms.comb.Nwds = 30;
parms.comb.Tw_opt = 'scale';
parms.comb.npol = 4;

parms.wind.pretime = -datwind(1);
parms.wind.prex = specwind(1);
parms.wind.postx = specwind(2);
parms.wind.taperx = 0.1;

parms.qc.minacor = 0.5;
parms.qc.maxphi = 5;

parms.inv.amp2phiwt = 2;
parms.inv.fmin = 0.15;
parms.inv.fmax = 0.5; % nominal - reset below according to fcross
parms.inv.corr_c_skip = true;
parms.inv.ifwt = true;

parms.inv.alpha = 0.25;

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

for ie = 216:402 % 44:norids % loop on orids
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
    elseif ~any([datinfo.xcor]==true)
        fprintf('Need to xcor arrival time for this phase and event\n'), continue
	end
    if ~isfield(datinfo,'dtstar') 
        fprintf('Need to calc specR for this phase and event\n'), continue
    elseif ~any([datinfo.dtstar]==true)
        fprintf('Need to calc specR for this phase and event\n'), continue
    end
    if isfield(datinfo,'comb')
        if any([datinfo.dtstar]==true)
            if ~overwrite
                yn = input('delta-tstar_comb already done - overwrite? [y/n] ','s'); 
                if ~strcmp(yn,'y'), fprintf('skipping...\n'), continue, end
            end
        end
    end
    
    if ifplotonly
        plot_ATTEN_TandF_domain_COMB( eqar )
        return
        continue
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

        % SHIFT TRACES USING Predicted ARRIVAL TIME 
        % shift so arrival is at time=0
        att = eqar(is).tt-eqar(is).pred_arrT; % shift to since absolute arrival
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
    
    % ONLY USE GOOD TRACES
    indgd = 1:size(eqar);
    indgd(mean(abs(all_dat0(:,indgd)))==0)     = []; % kill zero traces
    indgd(mean(abs(all_dat0(:,indgd)))<1e-12)     = []; % kill zero traces
    indgd(isnan(mean(abs(all_dat0(:,indgd))))) = []; % kill nan traces
    if ~isfield(eqar,'snr_wf'), continue; end % skip if snr not even calculated
    indgd([eqar(indgd).snr_wf]<snrmin)                  = []; % kill low snr traces
    if length(indgd) < 2, fprintf('NO GOOD TRACES/ARRIVALS, skip...\n'), continue, end

    %% Pick fmax from crossing freqs
    parms.inv.fmax = nanmean([eqar(indgd).fcross]');

    %% RUN THROUGH COMB
    [delta_tstar_comb,delta_T_comb,std_dtstar_comb,pairwise,fmids] = combspectra(all_dat0(:,indgd),resamprate,parms,0);

    %% ------------------ ALL-IN-ONE INVERSION + ALPHAS ------------------
    Amat = pairwise.As;
    phimat = pairwise.phis;
    wtmat = double(pairwise.inds).*pairwise.wts;
    % test_alphas = [0];
    parms.inv.amp2phiwt = 5;

%     %% use pairwise 1by1 then lsqr on each evt's dtstar & dT to get sta values
%     [ delta_tstar_pref,delta_T_pref,alpha_pref,alpha_misfits,dtstars,dTs ] ...
%         = invert_1by1_Aphis_4_STA_dtdtstar_alpha(Amat,phimat,fmids,test_alphas,wtmat,parms.inv.amp2phiwt );
%         
    %% use pairwise to make 1 global inversion for sta dtstar & dT
    [ delta_tstar_pref,delta_T_pref,alpha_pref,alpha_misfits,dtstars,dTs,A0s, Eamp, Ephi ] ...
        = invert_allin1_Aphis_4_STA_dtdtstar_alpha(Amat,phimat,fmids,test_alphas,wtmat,parms.inv.amp2phiwt );
    
    delta_tstar_use = delta_tstar_pref;
    delta_T_use = delta_T_pref;
    parms.inv.alpha = alpha_pref;
    
    %% ---------------------- STORE RESULTS -----------------------
    % STORE RESULTS
    fprintf('Recording results in arrival structure...')
    % prep eqar to receive new fields
    eqar(1).dtstar_comb = []; 
    eqar(1).dT_comb = []; 
%     eqar(1).stds_comb = []; 
    eqar(1).par_dtstar_comb = [];
                      
    eqar(indgd) =  dealto(eqar(indgd),'dtstar_comb',delta_tstar_use);
    eqar(indgd) =  dealto(eqar(indgd),'dT_comb',delta_T_use);
%     eqar(indgd) =  dealto(eqar(indgd),'stds_comb',std_dtstar_comb);
    eqar(indgd) =  dealto(eqar(indgd),'par_dtstar_comb',parms);

	indbd = setdiff(1:length(eqar),indgd);
    if ~isempty(eqar(indbd))
        eqar(indbd) =  dealto(eqar(indbd),'dtstar_comb',nan);
    end
    
	%% -------------------------- PLOTS ---------------------------
    if ifplot
    plot_ATTEN_TandF_domain_COMB( eqar )
    end % ifplot

%% -------------------------- SAVE ---------------------------
    if ifsave
    % SAVE
    save(arfile,'eqar')
    fprintf(' saved\n')
    % RECORD DTSTARCALC IN DATINFO 
    [datinfo.comb] = deal(false);
	[datinfo(indgd).comb] = deal(true);
    save(datinfofile,'datinfo')
    
    fprintf(' STA  CHAN  NEZ  resp  tilt  comp  xcor  specR  comb\n')
	for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f     %1.0f     %1.0f     %1.0f     %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).rmtilt,datinfo(is).rmcomp,datinfo(is).xcor,datinfo(is).dtstar,datinfo(is).comb); end
    end

    
    toc  
end % loop on orids

% results_PARSE