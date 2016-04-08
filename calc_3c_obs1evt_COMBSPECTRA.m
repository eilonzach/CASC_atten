% cycle through events and calculate spectral ratios for all stations
clear all
close all
cd /Users/zeilon/Documents/MATLAB/CASC_atten
addpath('matguts')

%% parameters
phase = 'S';
component = 'T'; %'Z', 'R', or 'T'
orid = 269;

samprate = 5 ; % new, common sample rate
filtfs = 1./[40 1]; % [flo fhi] = 1./[Tmax Tmin] in sec
taperx = 0.2;
datwind = [-160 165]; % window of data in eqar structure
specwind = [-5 30];
snrmin = 50;

unit = 'vel'; % 'disp'/'vel'/'acc'.

ifplot    = true;
ifsave    = false;

% combspec parms
parms.comb.Tmin = 1;
parms.comb.Tmax = 20;
parms.comb.Nwds = 30;
parms.comb.Tw_opt = 'scale';
parms.comb.npol = 4;

parms.wind.pretime = 200;
parms.wind.prex = -specwind(1);
parms.wind.postx = specwind(2);
parms.wind.taperx = 0.1;

parms.qc.minacor = 0.5;
parms.qc.maxphi = 5;

parms.inv.amp2phiwt = 2;
parms.inv.fmin = 0.15;
parms.inv.fmax = .5;
parms.inv.corr_c_skip = true;
parms.inv.ifwt = true;

%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% DATA DIRECTORY (top level)
datadir = '/Volumes/DATA_mini2/CASCADIA/DATA/'; % needs final slash

%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

%% GET EVENTS DATA
db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,elats,elons,edeps,evtimes,mags] = dbgetv(dbor,'orid','lat','lon','depth','time','ms');
norids = dbnrecs(dbor);
dbclose(db);

ie = find(orids==orid); 
evtime = evtimes(ie);
fprintf('\n Orid %.0f %s \n\n',orids(ie),epoch2str(evtime,'%Y-%m-%d %H:%M:%S'))

% name files and directories
evdir       = [num2str(orids(ie),'%03d'),'_',epoch2str(evtime,'%Y%m%d%H%M'),'/'];
datinfofile = [datadir,evdir,'_datinfo_',phase];
arfile      = [datadir,evdir,'_EQAR_',phase,'_',component];

% check files exist
if ~exist([datinfofile,'.mat'],'file'), fprintf('No station mat files for this event\n');return, end
if ~exist([arfile,'.mat'],'file'), fprintf('No arfile for this event + phase\n');return, end

% load files
load(datinfofile) % loads datinfo stucture
load(arfile)      % loads eqar structure

% options to skip
if isempty(datinfo), fprintf('No station mat files for this event\n'); return, end
if ~isfield(datinfo,'xcor') 
    fprintf('Need to xcor arrival time for this phase and event\n'), return
end
if ~any([datinfo.xcor]==true)
    fprintf('Need to xcor arrival time for this phase and event\n'), return
end

%% ONLY OBS STATIONS
db = dbopen([dbdir,dbnam],'r');
dbsi = dblookup_table(db,'site');
[stas,statypes] = dbgetv(dbsi,'sta','statype');
dbclose(db);
[~,inds,~] = intersect(stas,{eqar.sta},'stable');
statypes(inds);
%restrict to only obs
eqar = eqar(strcmp(statypes(inds),'OBS')); 
datinfo = datinfo(strcmp(statypes(inds),'OBS'));

%% ------------------ GRAB DATA IN RIGHT FORMAT ------------------
% prep data structures
nstas = length(datinfo);
all_dat0  = zeros(samprate*diff(datwind),nstas);

% LOOP ON STAS
for is = 1:nstas
    fprintf('Station %s... ',datinfo(is).sta)        % APPLY DATA QUALITY CONDITIONS
    if isempty(eqar(is).(['dat' component])),fprintf('no data for this component\n'),continue; end   
    if isempty(eqar(is).dT),fprintf('no xcor for this component\n'),continue; end
    if isnan(eqar(is).dT),fprintf('no xcor for this component\n'),continue; end
    if strcmp(eqar(is).sta,'M08C'), eqar(is).slon = -124.895400; end

    % SHIFT TRACES USING XCORRED ARRIVAL TIME 
    % shift so arrival is at time=0
    att = eqar(is).tt-eqar(is).abs_arrT; % shift to since MEASURED absolute arrival
    ja = (att >= datwind(1)) & (att < datwind(2)); % excerpt times according to datwind

    % GRAB DESIRED COMPONENT
    switch component
        case 'Z', all_dat0(:,is) = eqar(is).datZ(ja);
        case 'R', all_dat0(:,is) = eqar(is).datR(ja);
        case 'T', all_dat0(:,is) = eqar(is).datT(ja);
    end
    switch unit
        case 'disp'
            fprintf('Keeping in displacement... ');             
        case 'vel'
            fprintf('Diff to velocity... '); 
            all_dat0(:,is) = gradient(all_dat0(:,is),1./samprate);
            switch component
                case 'Z', eqar(is).datZ = gradient(eqar(is).datZ,1./eqar(is).samprate);
                case 'R', eqar(is).datR = gradient(eqar(is).datR,1./eqar(is).samprate);
                case 'T', eqar(is).datT = gradient(eqar(is).datT,1./eqar(is).samprate);
            end
        case 'acc'
            fprintf('Diff^2 to acceleration... '); 
            all_dat0(:,is) = gradient(gradient(all_dat0(:,is),1./samprate),1./samprate);
            switch component
                case 'Z', eqar(is).datZ = gradient(gradient(eqar(is).datZ,1./eqar(is).samprate),1./samprate);
                case 'R', eqar(is).datR = gradient(gradient(eqar(is).datR,1./eqar(is).samprate),1./samprate);
                case 'T', eqar(is).datT = gradient(gradient(eqar(is).datT,1./eqar(is).samprate),1./samprate);
            end
    end
    
    fprintf('got data\n')
end % loop on stas

% ONLY USE GOOD TRACES
indgd = 1:size(eqar);
indgd(mean(abs(all_dat0(:,indgd)))==0)     = []; % kill zero traces
indgd(mean(abs(all_dat0(:,indgd)))<1e-12)     = []; % kill zero traces
indgd(isnan(mean(abs(all_dat0(:,indgd))))) = []; % kill nan traces
indgd([eqar(indgd).snr_wf]<snrmin)                  = []; % kill low snr traces
if length(indgd) < 2, fprintf('NO GOOD TRACES/ARRIVALS, skip...\n'), return, end

%% MAKE STACK AND REFERENCE SPECTRUM
specss_ref = zeros(size(eqar(indgd(1)).specss));
for ig = 1:length(indgd)
    is = indgd(ig);
    specss_ref = specss_ref + eqar(is).specss;
end
specss_ref = specss_ref/length(indgd);

stak = sum(all_dat0(:,indgd),2)/length(indgd);
all_dat_do = [stak,all_dat0(:,indgd)];

%% Pick fmax from crossing freqs
parms.inv.fmax = nanmean([eqar(indgd).fcross]');

%% ------------------ COMB THE SPECTRA! ------------------
[delta_tstar_comb,delta_T_comb,std_dtstar_comb,pairwise,fmids] = combspectra_nofuss(all_dat_do,samprate,parms,false);
% [delta_tstar_comb_1,delta_T_comb_1,std_dtstar_comb_1,pairwise_1] = combspectra_cskip(all_dat_do,samprate,parms,false);
% [delta_tstar_comb_2,delta_T_comb_2,std_dtstar_comb_2,pairwise_2] = combspectra(all_dat_do,samprate,parms,false);

%% ------------------ ALL-IN-ONE INVERSION + ALPHAS ------------------
Amat = pairwise.As;
phimat = pairwise.phis;
wtmat = double(pairwise.inds).*pairwise.wts;
test_alphas = [0:0.05:0.9];
% test_alphas = [0];
parms.inv.amp2phiwt = 2;
 
[ delta_tstar_2,delta_T_2,alpha_pref,alpha_misfits,dtstars,dTs,A0s, Eamp, Ephi ] ...
    = allin1invert_Aphis_4_STA_dtdtstar_alpha(Amat,phimat,fmids,test_alphas,wtmat,parms.inv.amp2phiwt );

[delta_tstar_comb,delta_tstar_2]

return

%% ---------------------- STORE RESULTS -----------------------
% STORE RESULTS
fprintf('Recording results in arrival structure...')
% prep eqar to receive new fields
eqar(1).dtstar_comb = []; eqar(1).dT_comb = []; eqar(1).stds_comb = []; eqar(1).par_dtstar_comb = [];

eqar(indgd) =  dealto(eqar(indgd),'dtstar_comb',detrend(delta_tstar_comb(2:end),'constant'));
eqar(indgd) =  dealto(eqar(indgd),'stds_comb',std_dtstar_comb(2:end));
eqar(indgd) =  dealto(eqar(indgd),'dT_comb',detrend(delta_T_comb(2:end),'constant'));
eqar(indgd) =  dealto(eqar(indgd),'par_dtstar_comb',parms);

indbd = setdiff(1:length(eqar),indgd);
if ~isempty(eqar(indbd))
    eqar(indbd) =  dealto(eqar(indbd),'dtstar_comb',nan);
end

%% -------------------------- PLOTS ---------------------------
if ifplot
    fprintf('\fPlotting...\n')
% get noise-crossing freqs
fcross = [eqar(indgd).fcross]';
fmax = parms.inv.fmax;
%% plot simple spectra
figure(87), clf, set(gcf,'position',[440 0 1000 1400]), hold on
for ig = 1:length(indgd)
    is = indgd(ig);

    nyq = eqar(is).samprate/2;
    n4=length(find(eqar(is).frq<=0.8*nyq));
    hps = semilogx(eqar(is).frq(1:n4),log10(eqar(is).specs(1:n4)),'b','LineWidth',1.5); % << here using smoothed spectrum 
    hpn = semilogx(eqar(is).frq(1:n4),log10(eqar(is).specn(1:n4)),'k','LineWidth',1.5);
    set(hps,'color',colour_get(eqar(is).slon,max([eqar.slon]),min([eqar.slon])))
    % xlim(ax(1),[0 0.8*nyq]);
    xlabel('Hz'); 
    ylabel('amplitude, nm/Hz')
    % a=axis(ax(1));
    % % ylim_max=10^(floor(log10(max(specs(1:n4))))+1);
    % % ylim(ax(1),ylim_max.*[1.e-4 1])
    xlim([0 fmax]);
    ylim([-9,-2]);
end % loop on good stas

%% plot As, phis, spectral ratios, and the fits
figure(88), clf, set(gcf,'position',[440 0 1000 1400]), hold on
fmids = 1./logspace(log10(parms.comb.Tmin),log10(parms.comb.Tmax),parms.comb.Nwds)';

for ig = 1:length(indgd)
    is = indgd(ig);

    subplot(311), hold on
    plot(eqar(is).frq,log(eqar(is).specss./specss_ref),'-',...
            'color',colour_get(eqar(is).slon,max([eqar(indgd).slon]),min([eqar(indgd).slon]),flipud(jet)))
%     xlabel('\textbf{frequency (Hz)}','FontSize',18,'Interpreter','latex')
    set(gca,'Xscale','linear','xlim',[0.03 fmax+0.1],'ylim',[-3.5 3.5])
    ylabel('$\mathbf{\ln {(R)}}$','FontSize',18,'Interpreter','latex')
    title('specR','FontSize',18)

    subplot(312), hold on
    scatter(fmids,log(pairwise.As(ig,:)),110*pairwise.wts(ig,:),'o','Markeredgecolor','none','MarkerFaceColor',...
            colour_get(eqar(is).slon,max([eqar(indgd).slon]),min([eqar(indgd).slon]),flipud(jet)))
    plot(fmids,log(pairwise.As(ig,:)),...
            'color',colour_get(eqar(is).slon,max([eqar(indgd).slon]),min([eqar(indgd).slon]),flipud(jet)))
    set(gca,'Xscale','linear','xlim',[0.03 fmax+0.1],'ylim',[-3.5 3.5])
%     xlabel('\textbf{frequency (Hz)}','FontSize',18,'Interpreter','latex')
    ylabel('$\mathbf{\ln {(R)}}$','FontSize',18,'Interpreter','latex')
    title('COMB','FontSize',18)

	subplot(313), hold on
    scatter(fmids,pairwise.phis(ig,:),110*pairwise.wts(ig,:),'or','MarkerFaceColor',...
            colour_get(eqar(is).slon,max([eqar(indgd).slon]),min([eqar(indgd).slon]),flipud(jet)))
    plot(fmids,pairwise.phis(ig,:),':','color',...
            colour_get(eqar(is).slon,max([eqar(indgd).slon]),min([eqar(indgd).slon]),flipud(jet)))
    set(gca,'Xscale','linear','xlim',[0.03 fmax+0.1])
    xlabel('\textbf{frequency (Hz)}','FontSize',18,'Interpreter','latex'), 
    ylabel('$\mathbf{\Delta \phi}$ \textbf{(s)}','FontSize',18,'Interpreter','latex')
end % loop on good stas

figure(90)
scatter(delta_tstar_comb,-delta_T_comb,100./std_dtstar_comb)

eqar = plot_obs1evt_COMB( eqar);
% plot_ATTEN_TandF_domain_COMB( eqar )
end % ifplot

%% -------------------------- SAVE ---------------------------
if ifsave

%% spectra
figure(88),clf, set(gcf,'position',[400 400 580 540])
set(gca,'fontsize',14,'XTick',[0.05:0.05:0.25])

for is = 1:length(eqar)
    subplot(211), hold on
    scatter(fmids,log(pairwise.As(is,:)),110*pairwise.wts(is,:),'ok','Linewidth',0.5,'MarkerFaceColor',...
            colour_get(eqar(is).Xrdg,max([eqar.Xrdg]),min([eqar.Xrdg]),flipud(jet)))
    plot(fmids,log(pairwise.As(is,:)),...
            'color',colour_get(eqar(is).Xrdg,max([eqar.Xrdg]),min([eqar.Xrdg]),flipud(jet)))
    set(gca,'Xscale','linear','xlim',[0.03 fmax+0.1],'ylim',[-3.5 3.5],'fontsize',14,'box','on','linewidth',2)
    %     xlabel('\textbf{frequency (Hz)}','FontSize',18,'Interpreter','latex')
    ylabel('$\mathbf{\ln {(R)}}$','FontSize',18,'Interpreter','latex')

    subplot(212), hold on
    scatter(fmids,pairwise.phis(is,:),110*pairwise.wts(is,:),'ok','Linewidth',0.5,'MarkerFaceColor',...
            colour_get(eqar(is).Xrdg,max([eqar.Xrdg]),min([eqar.Xrdg]),flipud(jet)))
    plot(fmids,pairwise.phis(is,:),':','color',...
            colour_get(eqar(is).Xrdg,max([eqar.Xrdg]),min([eqar.Xrdg]),flipud(jet)))
    set(gca,'Xscale','linear','xlim',[0.03 fmax+0.1],'ylim',[-3.5 4.5],'fontsize',14,'box','on','linewidth',2)
    xlabel('\textbf{frequency (Hz)}','FontSize',18,'Interpreter','latex'), 
    ylabel('$\mathbf{\Delta \phi}$ \,\,\,\textbf{(s)}','FontSize',18,'Interpreter','latex')
end 

if ifsave
save2pdf(88,sprintf('1evt_OBS_spectra_%s_%s_%.0f_%s',...
    phase,component,orid,epoch2str(evtime,'%Y-%m-%d')),'figs')
end

%% waveforms
figure(2), 
set(gca,'fontsize',14)
if ifsave
    save2pdf(2,sprintf('1evt_OBS_waveforms_%s_%s_%.0f_%s',...
        phase,component,orid,epoch2str(evtime,'%Y-%m-%d')),'figs')
end

%% mapview
figure(32), 
% set(gcf,'position',[400 400 750 450])
% axis([lonlims latlims]+[-1 0 -1.4 1.5])
set(gca,'fontsize',14,'box','on','Linewidth',2)
title('$\Delta t^*$ recorded across JdF OBS stations','interpreter','latex','fontsize',18)

if ifsave
    save2pdf(32,sprintf('1evt_OBS_dtstar_map_%s_%s_%.0f_%s',...
        phase,component,orid,epoch2str(evtime,'%Y-%m-%d')),'figs')
end

%% section
figure(17), 
% title('Section of $\Delta t^*$ and $\delta T$ recorded across JdF ','interpreter','latex','fontsize',18)
if ifsave
    save2pdf(17,sprintf('1evt_OBS_section_%s_%s_%.0f_%s',...
        phase,component,orid,epoch2str(evtime,'%Y-%m-%d')),'figs')
end


end % ifsave

