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
clear all
close all
% close all
addpath('matguts')

% filter combs parms
Tmin = 1;
Tmax = 20;
% Twid = 2
Nwds = 30;

% window parms
pretime = 200;
prex = 10;
postx = 25;

minacor = 0.5;
maxphi = 5;
amp2phiwt = 2;

fmax = 1;

datdir = '/Volumes/DATA/CASCADIA/DATA/';
phase = 'S';
comp  = 'T';

% orid = 215; % seaz ~ 299  
orid = 230; % seaz ~ 299   
% orid = 238; % seaz ~ 230   
% orid = 269; % seaz ~ 250   

% sta1 = 'J53C'; % close to trench
% sta1 = 'J43C'; % close to trench
sta1 = 'J34C'; % close to trench

sta2 = 'J31C'; % close to ridge
% sta2 = 'J39C'; % close to ridge
% sta2 = 'J47C'; % close to ridge

%% ========== Load the two traces ==========
evt = dir([datdir,num2str(orid),'*']);
load([datdir,evt.name,'/_EQAR_',phase,'_',comp])
% load('synthetic_tstar/eg_eqar_456ST.mat')

ind1 = find(strcmp({eqar.sta},sta1)); %close to trench %  1  1  1  2
ind2 = find(strcmp({eqar.sta},sta2)); %close to ridge % 48 -3  5 30

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
fmax = mean([eqar([ind1,ind2]).fcross])

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

[ dat1 ] = filt_quick( dat1,1./30,10,dt);
[ dat2 ] = filt_quick( dat2,1./30,10,dt);

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
    
    figure(2), clf, set(gcf,'pos',[30 350 1200,400]), hold on
    plot(ttws,qdatwf2,'b','LineWidth',1)
    plot(ttws,qdatwf2s,'r','LineWidth',1)
    plot(ttws,qdatwf1,'k','LineWidth',2)
    plot(ttws,qdatwf2sa,'g','LineWidth',1)
    
    title(sprintf('Period = %.2f, A$_{obs}$ = %.2f, $\\phi_{obs}$ = %.2f, acor = %.2f',...
        Tmids(ii),A_f_obs,phi_f_obs,acor),'Fontsize',18,'interpreter','latex')
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
[ dtstar_e3,dT_e3,A0_e3 ] = invert_Aphi_4_dtdtstar( As(inds),phis(inds),fmids(inds), wts(inds),amp2phiwt);

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
parms.inv.fmin = .15;
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


%% plot for the poster/paper
As = pairwise_skip.As';
phis = pairwise_skip.phis';
wts = pairwise_skip.wts';
fmids = 1./logspace(log10(Tmin),log10(Tmax),Nwds)';
inds = find(fmids<=fmax & abs(phis)<maxphi & sqrt(wts)>minacor);
[ dtstar,dT,A0,misfit,res ] = invert_Aphi_4_dtdtstar( As(inds),phis(inds),fmids(inds), wts(inds),amp2phiwt);
frq = eqar(ind1).frq;

figure(30), clf, set(gcf,'pos',[600 600 600,800])

% ================ PLOT THE WAVEFORMS ================
subplot(5,1,1), hold on
plot(tt1,dat1./max(abs(dat1)),'k','LineWidth',2)
plot(tt2,dat2./max(abs(dat2)),'r','LineWidth',2)
plot(-prex*[1 1],[-1 1],'b--',postx*[1 1],[-1 1],'b--','linewidth',1.5)
text(-prex-18,-0.5,'\textbf{Station 1}','interpreter','latex','fontsize',12,'horizontalalignment','left','verticalalignment','bottom','color','k')
text(-prex-18,-0.8,'\textbf{Station 2}','interpreter','latex','fontsize',12,'horizontalalignment','left','verticalalignment','bottom','color','r')
text(postx+17,-1.5,'Time $\rightarrow$','interpreter','latex','fontsize',15,'horizontalalignment','left','verticalalignment','bottom')


% xlabel('Time, s','FontSize',22,'interpreter','latex')
set(gca,'xlim',[-prex-20 postx+20],'ylim',[-1.1 1.1],...%'XAxisLocation','top',...
    'Fontsize',12,'linewidth',2,'box','on','ytick',[])

% ================ PLOT THE AMPLITUDE SPECTRA ================
subplot(5,1,2:3), hold on
% plot(frq(2:end),log(spec2./spec1),'-ok')
plot(eqar(ind1).frq,log(eqar(ind2).specs./eqar(ind1).specs),'-xk','linewidth',1.5)
% plot(eqar(ind1).frq,log(eqar(ind2).specs./eqar(ind2).specn),'-ob')
% plot(eqar(ind1).frq,log(eqar(ind1).specs./eqar(ind1).specn),'-om')
plot(frq,log(A0) - pi*frq*dtstar,'g','Linewidth',2.5)

scatter(fmids,log(As),150*wts,'or','MarkerFaceColor','r')
scatter(fmids(inds),log(As(inds)),150*wts(inds),'ok','MarkerFaceColor','r','linewidth',1.5)

ha1 = arrow([eqar(ind1).fcross,1],[eqar(ind1).fcross,0],1.4,90,6,0.12,'FaceColor','m'); % fcross for sta 1
ha2 = arrow([eqar(ind2).fcross,1],[eqar(ind2).fcross,0],1.4,90,6,0.12,'FaceColor','g'); % fcross for sta 2
text(eqar(ind1).fcross,1.3,'$f^{max}_1$','interpreter','latex','fontsize',16,'horizontalalignment','center','verticalalignment','bottom')
text(eqar(ind2).fcross,1.3,'$f^{max}_2$','interpreter','latex','fontsize',16,'horizontalalignment','center','verticalalignment','bottom')
% plot(fmax*[1 1],[-1 1],'--b')

% xlabel('freq (Hz)','FontSize',16,'interpreter','latex')
ylabel('$\ln\,(R_{12})$','FontSize',22,'interpreter','latex')
set(gca,'fontsize',15,'Xscale','linear','xlim',[0.00 0.6],'ylim',[-3 3],'linewidth',2,'box','on')

% ================ PLOT THE PHASE SPECTRA ================
subplot(5,1,4:5), hold on
scatter(fmids,phis,150*wts,'or','MarkerFaceColor','r')
scatter(fmids(inds),phis(inds),150*wts(inds),'ok','MarkerFaceColor','r','linewidth',1.5)
plot(frq,dtstar*(1 - log(frq)./pi) + dT,'g','Linewidth',2.5)

% ha1 = arrow([eqar(ind1).fcross,1],[eqar(ind1).fcross,0],5,90,2,0.1,'FaceColor','m'); % fcross for sta 1
% ha2 = arrow([eqar(ind2).fcross,1],[eqar(ind2).fcross,0],5,90,2,0.1,'FaceColor','g'); % fcross for sta 2
% plot(fmax*[1 1],[-1 1],'--b')  
text(0.025,-1.6,['$\mathbf{\Delta t^* = ',num2str(dtstar,2),'}$'],'interpreter','latex','fontsize',20,'horizontalalignment','left','verticalalignment','bottom')


xlabel('frequency (Hz)','FontSize',22,'interpreter','latex')
ylabel('$\Delta \phi_{12}$ (s)','FontSize',22,'interpreter','latex')
set(gca,'fontsize',15,'Xscale','linear','xlim',[0.00 0.6],'ylim',[-2 3.5],'linewidth',2,'box','on')

% ================ SAVE FIGURE ================
ofile = sprintf('eg_comb_orid%.0f_%s%s_%s_v_%s',orid,phase,comp,sta1,sta2);
save2pdf(30,ofile,'figs');
