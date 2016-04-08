% Third script to test the method of measuring differential
% attenuation using combs of filters and measuring amplitude loss as well
% as phase lag within each narrow bandpass.
% This script uses real data and two different traces that have passed
% through different Q and V structure. 
% Second simple script to test the method of measuring differential
% attenuation using combs of filters and measuring amplitude loss as well
% as phase lag within each narrow bandpass.
% This script is attempts to do this with two traces that have passed
% through different Q and V material, to see how the velocity and Q
% differences trade off, and how the delta-tstar works.
clear all
% close all
addpath('matguts')

% filter combs parms
Tmin = 1;
Tmax = 20;
% Twid = 2
Nwds = 25;

% window parms
pretime = 200;
prex = 20;
postx = 30;

minacor = 0.5;
maxphi = 5;
amp2phiwt = 2;

fmax = 1;

%% ========== Load the two traces ==========
% load('synthetic_tstar/eg_eqar_456ST.mat')
% load('synthetic_tstar/eg_eqar_456ST.mat')
% load('/Volumes/DATA/CASCADIA/DATA/215_201308301625/_EQAR_S_T.mat') % seaz ~ 299  
load('/Volumes/DATA_mini2/CASCADIA/DATA/230_201310251710/_EQAR_S_T.mat') % seaz ~ 299   
% load('/Volumes/DATA/CASCADIA/DATA/238_201311230748/_EQAR_S_T.mat') % seaz ~ 230  
% load('/Volumes/DATA/CASCADIA/DATA/269_201404122014/_EQAR_S_T.mat')  %seaz ~ 250 

ind1 = find(strcmp({eqar.sta},'J53C')); %close to trench %  1  1  1  2
% ind1 = find(strcmp({eqar.sta},'J43C')); %close to trench %  1  1  1  2
ind1 = find(strcmp({eqar.sta},'J34C')); %close to trench %  1  1  1  2


ind2 = find(strcmp({eqar.sta},'J31C')); %close to ridge % 48 -3  5 30
% ind2 = find(strcmp({eqar.sta},'J39C')); %close to ridge % 48 -3  5 30
% ind2 = find(strcmp({eqar.sta},'J47C')); %close to ridge % 48 -3  5 30
ind2 = find(strcmp({eqar.sta},'J41C')); %close to ridge % 48 -3  5 30


dat1 = eqar(ind1).datT';
tt1 = eqar(ind1).tt'- eqar(ind1).pred_arrT;

if ~isstr(ind2)
    dat2 = eqar(ind2).datT';
    tt2 = eqar(ind2).tt'- eqar(ind2).pred_arrT;
end

samprate = eqar(ind1).samprate;
dt = 1./samprate;
fnq = samprate/2;
T = prex+postx;

%% MAKE STACK INSTEAD
if isstr(ind2)
indgd = find(~isnan([eqar.dT]));
tt2 = [-pretime:dt:pretime]'; tt2 = tt2(1:length(tt1));
for ig = 1:length(indgd)
    is = indgd(ig);
    att = eqar(is).tt-eqar(is).abs_arrT; % shift to since MEASURED absolute arrival
    all_dat0(:,is) = interp1(att,eqar(is).datT,tt2,'linear',0);
end % loop on stas
indgd = 1:size(eqar);
indgd(mean(abs(all_dat0(:,indgd)))<1e-6)     = []; % kill zero traces
indgd(isnan(mean(abs(all_dat0(:,indgd))))) = []; % kill nan traces
stak = sum(all_dat0(:,indgd),2)/length(indgd);
dat2 = stak;
end

%% CHANGE UNIT
% displacement to velocity
dat1_d = dat1;
dat2_d = dat2;
dat1_v = gradient(dat1_d,dt);
dat2_v = gradient(dat2_d,dt);
dat1_a = gradient(dat1_v,dt);
dat2_a = gradient(dat2_v,dt);
dat1 = dat1_a;
dat2 = dat2_a;
% 


%% filter 
[ dat1 ] = filt_quick( dat1,1./40,2,dt);
[ dat2 ] = filt_quick( dat2,1./40,2,dt);

%% Plot two traces
figure(1), clf, set(gcf,'pos',[100 550 1000 350]), hold on
plot(tt1,dat1./max(abs(dat1)),'k','LineWidth',1.5)
plot(tt2,dat2./max(abs(dat2)),'r','LineWidth',1.5)
plot(-prex*[1 1],max(abs(dat1))*[-1 1],'b--',postx*[1 1],max(abs(dat1))*[-1 1],'b--')


%% Make set of period windows for bandpass filter
Tmids = logspace(log10(Tmin),log10(Tmax),Nwds)';
Twdhs = 0.5*diff(logspace(log10(Tmin/2),log10(2*Tmax),Nwds+1)');
% Twdhs = 0.5*linspace(1,3,Nwds);
fmids = 1./Tmids;

As = zeros(Nwds,1);
phis = zeros(Nwds,1);
wts = zeros(Nwds,1);
for ii = 1:Nwds
    flo = 1./(Tmids(ii) + Twdhs(ii));
    fhi = 1./(Tmids(ii) - Twdhs(ii));
    fmid = 1./Tmids(ii);
    cp = struct('samprate',samprate,'pretime',pretime,'prex',prex,'postx',postx,...
                'taperx',0.1,'fhi',fhi,'flo',flo,'npoles',2,'norm',0);
    
	[ qdatwf1, qdatf1, qdatwc1, ~, ~, ttws, ~ ] = data_clean( dat1,cp );
    [ qdatwf2, qdatf2, qdatwc2, ~, ~, ~, ~ ] = data_clean( dat2,cp );

    % find observed phase shift
    [dcor, dcstd, dvcstd, acor]=xcortimes([qdatwf1,qdatwf2], dt, pretime, 10,0);
    phi_f_obs = diff(dcor);
    
    % make phase-corrected time series
    qdatwf2s = interp1(ttws-phi_f_obs,qdatwf2,ttws,'linear',0)';
    
    % calc. observed amplitude factor
    A_f_obs = (qdatwf2s'*qdatwf1)/(qdatwf1'*qdatwf1);
%     A_f_obs = 1./((datwf'*qdatwfs)/(datwf'*datwf)); %<< could do it other way too
    
    % make amplitude-corrected time series
    qdatwf2sa = qdatwf2s./A_f_obs;
    
    acor = xcorr(qdatwf1,qdatwf2sa,0)^2./(xcorr(qdatwf1,qdatwf1,0)*xcorr(qdatwf2sa,qdatwf2sa,0));
    
%     figure(2), clf, set(gcf,'pos',[30 350 1200,400]), hold on
%     plot(ttws,qdatwf2,'b','LineWidth',1)
%     plot(ttws,qdatwf2s,'r','LineWidth',1)
%     plot(ttws,qdatwf1,'k','LineWidth',2)
%     plot(ttws,qdatwf2sa,'g','LineWidth',1)
%     
%     title(sprintf('Period = %.2f, A$_{obs}$ = %.2f, $\\phi_{obs}$ = %.2f, acor = %.2f',...
%         Tmids(ii),A_f_obs,phi_f_obs,acor),'Fontsize',18,'interpreter','latex')
%     xlim([-30 30])

%       pause

    As(ii) = A_f_obs;
    phis(ii) = phi_f_obs;
    wts(ii) = acor.^2;
end



%% use vectors of amplitude factors and phase shifts to estimate dtstar and dT
% in theory, A(f)   = exp(pi*f*tstar)
%  so a plot of ln(A) against freq should have slope pi*dtstar
%  or a plot of ln(A) against 1/t should have slope pi*dtstar
% in theory, phi(f) = 1./pi * tstar * ln(fnq/f)
%  so a plot of phi against ln(freq) should have slope -dtstar./pi and
%  intercept that related to fnq and dt
%  or a plot of phi against ln(t) should have slope dtstar./pi
figure(4), clf, set(gcf,'pos',[600 10 500,700])
subplot(211), hold on
scatter(fmids,log(As),110*wts,'or','MarkerFaceColor','r'), set(gca,'Xscale','log')
xlabel('log(freq)','FontSize',18), ylabel('log(Amp)','FontSize',18)
subplot(212), hold on
scatter(fmids,phis,110*wts,'or','MarkerFaceColor','r'), set(gca,'Xscale','log')
xlabel('log(freq)','FontSize',18), ylabel('phase-shift (s)','FontSize',18)
 
%% QC
inds = find(fmids<=fmax & abs(phis)<maxphi & sqrt(wts)>minacor);

% estimates of dtstar from the data
fo1 = fit(fmids(inds),log(As(inds)),'poly1','weight',wts(inds));
dtstar_e1 = -fo1.p1./pi;
fo2 = fit(log(fmids(inds)),phis(inds),'poly1','weight',wts(inds));
dtstar_e2 = -fo2.p1*pi;
% estimates of dT from the data
dT_anel = fo2.p2 + fo2.p1*pi; % using phase lags, taking into account anelasticity from phase

figure(5);plot([qdatwc1,qdatwc2])
dT_xcor = diff(xcortimes([qdatwc1,qdatwc2], dt, prex, maxphi,1)); % ignoring anelasticity

% estimates from simultaneous inversion of amp and phase data
a_tests = [0:0.05:0.9]';
a_misfits = zeros(length(a_tests),1);
for ia = 1:length(a_tests)
[ ~,~,~,a_misfits(ia),~ ] = invert_Aphi_4_dtdtstar( As(inds),phis(inds),fmids(inds), wts(inds),1,a_tests(ia));
end
figure(7), clf, hold on
plot(a_tests,a_misfits)
fs = fit(a_tests,a_misfits,'smoothingspline');
plot(fs), xlabel('alpha'), ylabel('misfit')
[ dtstar_e3,dT_e3,A0_e3] = invert_Aphi_4_dtdtstar( As(inds),phis(inds),fmids(inds), wts(inds),3,a_tests(mindex(a_misfits)));

%% estimate dtstar from spectral ratio
wlen = length(qdatwc1); % window length, accounting for taper+padding
ntap=2;
nft=2^nextpow2(wlen);
[spec1,frq]=pmtm(qdatwc1,ntap,nft,samprate);
[spec2,~]=pmtm(qdatwc2,ntap,nft,samprate);
spec1 = spec1(2:length(spec1)).^0.5;
spec2 = spec2(2:length(spec2)).^0.5;
figure(4); clf
dtstar_lnR =  diff(xspecratio( [spec1,spec2],frq,fmax,0.01,0,1 ));


%% plot estimated vals back on, to test
figure(4)
subplot(211), hold on
plot(frq(2:end),log(spec2./spec1),'-ok')
plot(fmids,fo1.p2 - pi*fmids*dtstar_e1,'b','Linewidth',2)
plot(fmids,-pi*fmids*dtstar_e2,'m','Linewidth',2)
plot(fmids,log(A0_e3) - pi*fmids*dtstar_e3,'g','Linewidth',2)
plot(fmax*[1 1],[-1 1],'--b','Linewidth',2)
scatter(fmids,log(As),110*wts,'or','MarkerFaceColor','r')% data again
% set(gca,'Xscale','log','xlim',[0.03 1.])
set(gca,'Xscale','linear','xlim',[0.03 0.6])
xlabel('freq','FontSize',18),ylabel('log(Amp)','FontSize',18)

subplot(212), hold on
plot(fmids,dtstar_e2*(1 - log(fmids)/pi) + dT_anel,'m'), set(gca,'Xscale','log')
plot(fmids,dtstar_e1*(1 - log(fmids)/pi),'b'), set(gca,'Xscale','log')
plot(fmids,dtstar_e3*(1 - log(fmids)/pi) + dT_e3,'g'), set(gca,'Xscale','log')
plot(fmax*[1 1],[-1 1],'--b','Linewidth',2)
scatter(fmids,phis,110*wts,'or','MarkerFaceColor','r') % data again
set(gca,'Xscale','log','xlim',[0.03 fmax],'ylim',maxphi*[-1 1])
xlabel('log(freq)','FontSize',18), ylabel('phase-shift (s)','FontSize',18)


%% summarise
fprintf('\nPositive dt* means 2 (red) more attenuated than 1 (black)\n')
fprintf('dt* est1 = %.3f from amplitude (blue)\n',dtstar_e1)
fprintf('dt* est2 = %.3f from phase (magenta)\n',dtstar_e2)
fprintf('dt* est3 = %.3f from both\n',dtstar_e3)
fprintf('dt* est3 = %.3f from spectral ratio\n\n',dtstar_lnR)

fprintf('\nPositive dT means 2 (red) slower than 1 (black)\n')
fprintf('dT est1\t= %.3f with anelasticity (phi)\n',dT_anel)
fprintf('dT est2\t= %.3f with anelasticity both\n',dT_e3)
fprintf('dT est2\t= %.3f xcorr assuming elastic\n',dT_xcor)

return

%% Test combspectra function
parms.comb.Tmin = Tmin; % fnq/2
parms.comb.Tmax = Tmax;
parms.comb.Nwds = Nwds;
parms.comb.Tw_opt = 'scale';
parms.comb.npol = 4;

parms.wind.pretime = pretime;
parms.wind.prex = prex;
parms.wind.postx = postx;
parms.wind.taperx = 0.1;

parms.qc.minacor = minacor;
parms.qc.maxphi = maxphi;

parms.inv.amp2phiwt = amp2phiwt;
parms.inv.fmin = .1;
parms.inv.fmax = fmax;
parms.inv.ifwt = true;
parms.inv.corr_c_skip = true;

ifplot = false;
% 
% tic
% [delta_tstar,delta_T,std_dtstar,pairwise] = combspectra([dat1,dat2],samprate,parms,1);
% fprintf('combspectra dtstar est\t= %.3f \n',diff(delta_tstar))
% fprintf('combspectra dT est\t= %.3f \n',diff(delta_T))
% toc
tic
[delta_tstar_skip,delta_T_skip,std_dtstar_skip,pairwise_skip] = combspectra_nofuss([dat1,dat2],samprate,parms,1);
fprintf('BETTER? combspectra dtstar est\t= %.3f \n',diff(delta_tstar_skip))
fprintf('BETTER? combspectra dT est\t= %.3f \n',diff(delta_T_skip))
toc


% %% Cross spectra
% 
% [Cxy,F] = cpsd(dat,qdat,50,0,[],samprate);
% % [Cxy,F] = cpsd(dat,dat,50,0,[],samprate);
% figure(99), clf
% plot(F,-angle(Cxy)/pi), xlim([0,2])
% figure(98), clf
% plot(F,abs(Cxy)), xlim([0,2])
% % plot(F,abs(angle(Cxy)))
% return
% 
% %% Plot filtered traces with and without extra attenuation
% for iff = 1:Nfs
%     filtfs = filtcomb(iff,:);
% 
%     cp = struct('samprate',samprate,'pretime',200,'prex',20,'postx',30,...
%                 'taperx',0.1,'fhi',filtfs(2),'flo',filtfs(1),'npoles',2,'norm',0);
%     [ datwf, datf, datwc, datc, ~, ttws, tts ] = data_clean( dat,cp );
%     [ qdatwf, qdatf, qdatwc, adatc, ~, ~, ~ ] = data_clean( qdat,cp );
%     
%     datr = datr+datf;
%     qdatr = qdatr+qdatf;
% 
%     figure(2), hold on
%     plot(ttws,datwc,'r','LineWidth',2)
%     plot(ttws,qdatwc,'--r','LineWidth',2)
%     plot(ttws,datwf,'color',colour_get(iff,Nfs,1))
%     plot(ttws,qdatwf,'--','color',colour_get(iff,Nfs,1))
%     plot(tts,datr,'k')
%     plot(tts,qdatr,'.k')
%     xlim([min(ttws),max(ttws)])
% 
%     
%     amps(iff) = norm(datwf)/delf;
%     qamps(iff) = norm(qdatwf)/delf;
% 
%     figure(3), hold on
%     plot(mean(filtfs),amps(iff)/qamps(iff)/mean(filtfs),'ob')
% end
% return
% 
% %% Plot filtered traces and their sum to show building the full signal
% for iff = 1:Nfs
%     filtfs = filtcomb(iff,:);
% 
%     cp = struct('samprate',samprate,'pretime',200,'prex',20,'postx',30,...
%                 'taperx',0.1,'fhi',filtfs(2),'flo',filtfs(1),'npoles',2,'norm',0);
%     [ datwf, datf, datwc, datc, ~, ttws, tts ] = data_clean( dat,cp );
%     
%     datr = datr+datf;
% 
%     figure(2), hold on
%     plot(ttws,datwc,'r','LineWidth',2)
%     plot(ttws,datwf,'color',colour_get(iff,Nfs,1))
%     plot(tts,datr,'k','color',colour_get(iff,Nfs,1))
%     
%     xlim([min(ttws),max(ttws)])
% 
% end
% return
% 
% 
% %% Plotting amplitudes of the filtered traces
% amps = zeros(Nfs,1);
% for iff = 1:Nfs
%     
%     filtfs = filtcomb(iff,:);
% 
%     cp = struct('samprate',samprate,'pretime',10,'prex',10,'postx',10,...
%                 'taperx',0.0,'fhi',filtfs(2),'flo',filtfs(1),'npoles',2,'norm',0);
%     [ datwf, datf, datwc, datc, ~, ttws, tts ] = data_clean( qdat,cp );
%     
%     
%     amps(iff) = norm(datwf)/delf;
% 
%     figure(2), hold on
%     plot(mean(filtfs),amps(iff),'ob')
%     
% end
% plot(ff,abs(DAT.*Dwt),'xr')
% 
% return
% 
% % load('synthetic_tstar/eg_eqar_456ST.mat')
% % dat0 = eqar(1).datT;
% return
% 
