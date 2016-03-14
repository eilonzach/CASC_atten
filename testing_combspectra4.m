% Third script to test the method of measuring differential
% attenuation using combs of filters and measuring amplitude loss as well
% as phase lag within each narrow bandpass.
% This script uses real data and three different traces that have passed
% through different Q and V structure. 
% Second simple script to test the method of measuring differential
% attenuation using combs of filters and measuring amplitude loss as well
% as phase lag within each narrow bandpass.
% This script is attempts to do this with two traces that have passed
% through different Q and V material, to see how the velocity and Q
% differences trade off, and how the delta-tstar works.
% clear all
% close all
cd /Users/zeilon/Documents/MATLAB/CASC_atten;
addpath('matguts')

parms.comb.Tmin = 1;
parms.comb.Tmax = 20;
parms.comb.Nwds = 30;
parms.comb.Tw_opt = 'scale';
parms.comb.npol = 2;

parms.wind.pretime = 200;
parms.wind.prex = 15;
parms.wind.postx = 25;
parms.wind.taperx = 0.1;

parms.qc.minacor = 0.5;
parms.qc.maxphi = 5;

parms.inv.amp2phiwt = 3;
parms.inv.fmin = 0.15;
parms.inv.fmax = 1
parms.inv.ifwt = true;

ifplot = true;


%% ========== Load the two traces ==========
% load('synthetic_tstar/eg_eqar_456ST.mat')
load('synthetic_tstar/eg_eqar_456ST.mat')
load('/Volumes/DATA/CASCADIA/DATA/269_201404122014/_EQAR_S_T.mat')  

dat = [];
stas = {};
for is = 1:40
    if isempty(eqar(is).datT), continue, end
    dat = [dat,eqar(is).datT'];
    stas = {stas{:},eqar(is).sta};
end

% correct bad chans for 456ST

tt = eqar(1).tt-eqar(1).pred_arrT;
samprate = eqar(1).samprate;

%% use filter comb to estimate delta_tstar and delta_T
[delta_tstar,delta_T,std_dtstar,pairwise] = combspectra_cskip(dat,samprate,parms,0)

%% use spectral ratios to estimate 
taperx = 0.2;
cp = struct('samprate',samprate,'pretime',parms.wind.pretime,'prex',parms.wind.prex,'postx',parms.wind.postx,...
            'taperx',taperx,'fhi',1./parms.comb.Tmin,'flo',1./parms.comb.Tmax,'npoles',2,'norm',0);
[ all_datwf,all_datf,all_datwc,all_datc,~,ttws,tts ] = data_clean( dat,cp );
wlen = (parms.wind.prex+parms.wind.postx)*samprate*(1+4*taperx); % window length, accounting for taper+padding
dt = 1./samprate;
ntap=2;
nft=2^nextpow2(wlen);
nyq=0.5*samprate;
specss =[];
hw = waitbar(0,'PROGRESS THROUGH CALCULATION OF SPECTRA');
for is = 1:size(dat,2)        
    [specs,frq]=pmtm(all_datwf(:,is),ntap,nft,samprate);
    frq=frq(2:length(frq));
    specs=specs(2:length(specs));

    %convert power to amplitude ALREADY IN DISPLACEMENT
    specss(:,is)=specs.^0.5;
%     eqar(is).specss = moving_average(eqar(is).specs,mavwind);

    waitbar(is/size(dat,2),hw)
end
delete(hw)

% CALC DELTA-TSTAR.
[ delta_tstar_specratio,cov_dtstar_specratio,std_dtstar_specratio ] = xspecratio( specss,frq,parms.inv.fmax,0,1,ifplot );

%% use cross-correlation to estimate anharmonic delta_T
% ind = tt>=-parms.wind.prex & tt<parms.wind.postx;
% [dcor, dcstd, dvcstd, acor]=xcortimes(dat(ind,:), 1./samprate, parms.wind.pretime, 5,1);
[dcor, dcstd, dvcstd, acor]=xcortimes(all_datwf, 1./samprate, parms.wind.pretime, 5,1);

figure(47), clf, hold on
scatter(delta_tstar,delta_T,std_dtstar*70,'o'), xlabel('delta tstar','Fontsize',16),ylabel('delta T anharm','Fontsize',16)
text(delta_tstar,delta_T,stas)

figure(48),clf, hold on
scatter(dcor,delta_T,'o'), xlabel('dcor','Fontsize',16),ylabel('delta T anel','Fontsize',16)

figure(49),clf, hold on
scatter(dcor,delta_tstar,'o'), xlabel('dcor','Fontsize',16),ylabel('delta tstar','Fontsize',16)



figure(50),clf, hold on
plot(dcor,delta_tstar_specratio,'o'), xlabel('dcor','Fontsize',16),ylabel('delta tstar specratio','Fontsize',16)

figure(51),clf
plot(delta_tstar,delta_tstar_specratio,'o'), xlabel('delta tstar comb','Fontsize',16),ylabel('delta tstar specratio','Fontsize',16)


% tt2 = eqar(ind2).tt'- eqar(ind2).pred_arrT;
% 
% dt = 1./samprate;
% fnq = samprate/2;
% T = prex+postx;
% 
% %% Plot two traces
% figure(1), clf, set(gcf,'pos',[100 550 1000 350]), hold on
% plot(tt1,dat1,'k','LineWidth',1.5)
% plot(tt2,dat2,'r','LineWidth',1.5)
% plot(-prex*[1 1],max(abs(dat1))*[-1 1],'b--',postx*[1 1],max(abs(dat1))*[-1 1],'b--')
% 
% %% Make set of period windows for bandpass filter
% Tmids = logspace(log10(Tmin),log10(Tmax),Nwds)';
% Twdhs = 0.5*diff(logspace(log10(Tmin/2),log10(2*Tmax),Nwds+1)');
% % Twdhs = 0.5*Tmids;
% fmids = 1./Tmids;
% fprintf('SHOULD PLOT FILTERS IN F-SPACE\n')
% 
% As = zeros(Nwds,1);
% phis = zeros(Nwds,1);
% wts = zeros(Nwds,1);
% for ii = 1:Nwds
%     flo = 1./(Tmids(ii) + Twdhs(ii));
%     fhi = 1./(Tmids(ii) - Twdhs(ii));
%     fmid = 1./Tmids(ii);
%     cp = struct('samprate',samprate,'pretime',pretime,'prex',prex,'postx',postx,...
%                 'taperx',0.1,'fhi',fhi,'flo',flo,'npoles',2,'norm',0);
%     
% 	[ qdatwf1, qdatf1, qdatwc1, ~, ~, ttws, ~ ] = data_clean( dat1,cp );
%     [ qdatwf2, qdatf2, qdatwc2, ~, ~, ~, ~ ] = data_clean( dat2,cp );
% 
%     % find observed phase shift
%     [dcor, dcstd, dvcstd, acor]=xcortimes([qdatwf1,qdatwf2], dt, pretime, 10,0);
%     phi_f_obs = diff(dcor);
%     
%     % make phase-corrected time series
%     qdatwf2s = interp1(ttws-phi_f_obs,qdatwf2,ttws,'linear',0)';
%     
%     % calc. observed amplitude factor
%     A_f_obs = (qdatwf2s'*qdatwf1)/(qdatwf1'*qdatwf1);
% %     A_f_obs = 1./((datwf'*qdatwfs)/(datwf'*datwf)); %<< could do it other way too
%     
%     % make amplitude-corrected time series
%     qdatwf2sa = qdatwf2s./A_f_obs;
%     
%     acor = xcorr(qdatwf1,qdatwf2sa,0)^2./(xcorr(qdatwf1,qdatwf1,0)*xcorr(qdatwf2sa,qdatwf2sa,0));
%     
%     figure(2), clf, set(gcf,'pos',[30 350 1200,400]), hold on
%     plot(ttws,qdatwf2,'b','LineWidth',1)
%     plot(ttws,qdatwf2s,'r','LineWidth',1)
%     plot(ttws,qdatwf1,'k','LineWidth',2)
%     plot(ttws,qdatwf2sa,'g','LineWidth',1)
%     
%     title(sprintf('Period = %.2f, A$_{obs}$ = %.2f, $\\phi_{obs}$ = %.2f, acor = %.2f',...
%         Tmids(ii),A_f_obs,phi_f_obs,acor),'Fontsize',18,'interpreter','latex')
% %     xlim([-30 30])
% 
%      pause
% 
%     As(ii) = A_f_obs;
%     phis(ii) = phi_f_obs;
%     wts(ii) = acor.^2;
% end
% 
% %% use vectors of amplitude factors and phase shifts to estimate dtstar and dT
% % in theory, A(f)   = exp(pi*f*tstar)
% %  so a plot of ln(A) against freq should have slope pi*dtstar
% %  or a plot of ln(A) against 1/t should have slope pi*dtstar
% % in theory, phi(f) = 1./pi * tstar * ln(fnq/f)
% %  so a plot of phi against ln(freq) should have slope -dtstar./pi and
% %  intercept that related to fnq and dt
% %  or a plot of phi against ln(t) should have slope dtstar./pi
% figure(3), clf, set(gcf,'pos',[600 10 500,700])
% subplot(211), hold on
% scatter(fmids,log(As),70*wts,'or','MarkerFaceColor','r'), set(gca,'Xscale','log')
% xlabel('log(freq)','FontSize',18), ylabel('log(Amp)','FontSize',18)
% subplot(212), hold on
% scatter(fmids,phis,70*wts,'or','MarkerFaceColor','r'), set(gca,'Xscale','log')
% xlabel('log(freq)','FontSize',18), ylabel('phase-shift (s)','FontSize',18)
% 
% fmax = input('What is fmax? ');
%  
% %% QC
% inds = fmids<=fmax & abs(phis)<maxphi & sqrt(wts)>minacor;
% 
% 
% % estimates of dtstar from the data
% fo1 = fit(fmids(inds),log(As(inds)),'poly1','weight',wts(inds));
% dtstar_e1 = -fo1.p1./pi;
% fo2 = fit(log(fmids(inds)),phis(inds),'poly1','weight',wts(inds));
% dtstar_e2 = -fo2.p1*pi;
% % estimates of dT from the data
% dT_anel = fo2.p2 + fo2.p1*log(fnq); % using phase lags, taking into account anelasticity from phase
% dT_xcor = diff(xcortimes([qdatwc1,qdatwc2], dt, 50, 10,0)); % ignoring anelasticity
% 
% % estimates from simultaneous inversion of amp and phase data
% [ dtstar_e3,dT_e3,A0_e3 ] = invert_Aphi_4_dtdtstar( As(inds),phis(inds),fmids(inds),fnq, wts(inds),amp2phiwt);
% 
% 
% % plot estimated vals back on, to test
% figure(3)
% subplot(211), hold on
% plot(fmids,fo1.p2 - pi*fmids*dtstar_e1,'b'), set(gca,'Xscale','log')
% plot(fmids,-pi*fmids*dtstar_e2,'m'), set(gca,'Xscale','log')
% plot(fmids,log(A0_e3) - pi*fmids*dtstar_e3,'g'), set(gca,'Xscale','log')
% 
% plot(fmax*[1 1],[-1 1],'--b')
% subplot(212), hold on
% plot(fmids,(log(fnq) - log(fmids))*dtstar_e2./pi + dT_anel,'m'), set(gca,'Xscale','log')
% plot(fmids,(log(fnq) - log(fmids))*dtstar_e1./pi,'b'), set(gca,'Xscale','log')
% plot(fmids,(log(fnq) - log(fmids))*dtstar_e3./pi + dT_e3,'g'), set(gca,'Xscale','log')
% plot(fmax*[1 1],[-1 1],'--b')
% 
% %% estimate dtstar from spectral ratio
% wlen = length(qdatwc1); % window length, accounting for taper+padding
% ntap=2;
% nft=2^nextpow2(wlen);
% [spec1,frq]=pmtm(qdatwc1,ntap,nft,samprate);
% [spec2,~]=pmtm(qdatwc2,ntap,nft,samprate);
% spec1 = spec1(2:length(spec1)).^0.5;
% spec2 = spec2(2:length(spec2)).^0.5;
% figure(4); clf
% dtstar_lnR =  diff(xspecratio( [spec1,spec2],frq,fmax,0.01,0,1 ));
% 
% %% summarise
% fprintf('\nPositive dt* means 2 more attenuated than 1\n')
% fprintf('dt* est1 = %.3f from amplitude (blue)\n',dtstar_e1)
% fprintf('dt* est2 = %.3f from phase (magenta)\n',dtstar_e2)
% fprintf('dt* est3 = %.3f from both\n',dtstar_e3)
% fprintf('dt* est3 = %.3f from spectral ratio\n\n',dtstar_lnR)
% 
% fprintf('\nPositive dT means 2 slower than 1\n')
% fprintf('dT est1\t= %.3f with anelasticity (amp)\n',dT_anel)
% fprintf('dT est2\t= %.3f with anelasticity both\n',dT_e3)
% fprintf('dT est2\t= %.3f xcorr assuming elastic\n',dT_xcor)
% return
% 
% 