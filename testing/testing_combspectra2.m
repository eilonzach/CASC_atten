% Second simple script to test the method of measuring differential
% attenuation using combs of filters and measuring amplitude loss as well
% as phase lag within each narrow bandpass.
% This script is attempts to do this with two traces that have passed
% through different Q and V material, to see how the velocity and Q
% differences trade off, and how the delta-tstar works.
% clear all
close all
addpath('matguts')
% parms
samprate = 100;
%% trace 1 Q & V - faster c0 and less attenuated
Q0_1 = 150;
c0_1 = 3.95e3; % reference velocity in km/s
%% trace 2 Q & V - slower c0 and more attenuated
Q0_2 = 50;
c0_2 = 3.95e3; % reference velocity in km/s

L = 200e3;

alpha = 0;

% calc. true values
dT = L*(1./c0_2 - 1./c0_1);
dtstar = L*(1./c0_2./Q0_2 - 1./c0_1./Q0_1);
fprintf('Theoretical dT=%.2f and dtstar=%.2f\n',dT,dtstar)

%% ========== Make the two traces ==========

%% Make simple pulse
dt = 1/samprate;
tt = [-50:dt:50]';
T = tt(end)-tt(1);
N = T/dt; tt = tt(1:N);
fnq = 0.5/dt;
dat = synthtrace(T,2,1,dt,'gauss');

%% Delay pulses
dT1 = L*(1./c0_1 - 1./mean([c0_1,c0_2]));
dat1 = interp1(tt+dT1,dat,tt,'linear',0)';
dT2 = L*(1./c0_2 - 1./mean([c0_1,c0_2]));
dat2 = interp1(tt+dT2,dat,tt,'linear',0)';

%% Attenuate pulses
% take fft of pulse
[DAT_1,~ ] = fft_ze(dat1,dt);
[DAT_2,ff] = fft_ze(dat2,dt);
w = 2*pi*ff;
% Make attenuation operators
[ Dwt_1 ] = attenuation_operator( Q0_1,c0_1,L,w,alpha);
[ Dwt_2 ] = attenuation_operator( Q0_2,c0_2,L,w,alpha);
% apply attenuation operator and ifft
qdat1 = abs(ifft(DAT_1.*Dwt_1));
qdat2 = abs(ifft(DAT_2.*Dwt_2));


%% Plot original and attenuated pulses
figure(1), clf, hold on
plot(tt,dat,'k','LineWidth',1.5)
plot(tt,qdat1,'b','LineWidth',1.5)
plot(tt,qdat2,'r','LineWidth',1.5)
title('Original (black), trace1 (blue) and trace2 (red)')

%% Make set of period windows for bandpass filter
Tmax = T;
Tmin = 2/fnq;
Nwds = 20;
Tmids = logspace(log10(2/fnq),log10(T),Nwds)';
Twdhs = 0.5*diff(logspace(log10(1/fnq),log10(T/2),Nwds+1)');
fmids = 1./Tmids;
As = zeros(Nwds,1);
phis = zeros(Nwds,1);

for ii = 1:Nwds
    flo = 1./(Tmids(ii) + Twdhs(ii));
    fhi = 1./(Tmids(ii) - Twdhs(ii));
    fmid = 1./Tmids(ii);
    cp = struct('samprate',samprate,'pretime',50,'prex',50,'postx',50,...
                'taperx',0.1,'fhi',fhi,'flo',flo,'npoles',2,'norm',0);
    [ qdatwf1, qdatf1, qdatwc1, ~, ~, ttws, ~ ] = data_clean( qdat1,cp );
    [ qdatwf2, qdatf2, qdatwc2, ~, ~, ~, ~ ] = data_clean( qdat2,cp );

    % find observed phase shift
    [dcor, dcstd, dvcstd, acor]=xcortimes([qdatwf1,qdatwf2], dt, 50, 10,0);
    phi_f_obs = diff(dcor);
    % calc. predicted phase shift
     phi_f_pred = dtstar*log(fnq/fmid)./pi + dT;

    
    % make phase-corrected time series
    qdatwf2s = interp1(ttws-phi_f_obs,qdatwf2,ttws,'linear',0)';
    % calc. observed amplitude factor
    A_f_obs = (qdatwf2s'*qdatwf1)/(qdatwf1'*qdatwf1);
%     A_f_obs = 1./((datwf'*qdatwfs)/(datwf'*datwf)); %<< could do it other way too
    % calc. predicted amplitude factor from tstar
    A_f_pred = exp(-pi*fmid*dtstar);
    
    % make amplitude-corrected time series
    qdatwf2sa = qdatwf2s./A_f_obs;
    
    figure(2), clf, set(gcf,'pos',[30 350 1200,400]), hold on
    plot(ttws,qdatwf1,'k','LineWidth',2)
    plot(ttws,qdatwf2,'b','LineWidth',1.5)
    plot(ttws,qdatwf2s,'r','LineWidth',1.5)
    plot(ttws,qdatwf2sa,'g','LineWidth',1)
    
    title(sprintf('Period = %.2f, A$_{pred}$ = %.2f, A$_{obs}$ = %.2f, $\\phi_{pred}$ = %.2f, $\\phi_{obs}$ = %.2f',...
        Tmids(ii),A_f_pred,A_f_obs,phi_f_pred,phi_f_obs),'Fontsize',18,'interpreter','latex')
    xlim([-20 20])

    As(ii) = A_f_obs;
    phis(ii) = phi_f_obs;
%     pause
end
%% use vectors of amplitude factors and phase shifts to estimate dtstar and dT
% in theory, A(f)   = exp(pi*f*tstar)
%  so a plot of ln(A) against freq should have slope pi*dtstar
%  or a plot of ln(A) against 1/t should have slope pi*dtstar
% in theory, phi(f) = 1./pi * tstar * ln(fnq/f)
%  so a plot of phi against ln(freq) should have slope -dtstar./pi and
%  intercept that related to fnq and dt
%  or a plot of phi against ln(t) should have slope dtstar./pi
ff = logspace(log10(fmids(end)),log10(fnq),40);
figure(3), clf, set(gcf,'pos',[600 10 500,700])
subplot(211), hold on
plot(fmids,log(As),'or')
plot(ff,-pi*ff*dtstar,'b'), set(gca,'Xscale','log')
xlabel('log(freq)','FontSize',18), ylabel('log(Amp)','FontSize',18)
subplot(212), hold on
plot(fmids,phis,'or')
plot(ff,dtstar - log(ff)*dtstar./pi + dT,'b'), set(gca,'Xscale','log')
xlabel('log(freq)','FontSize',18), ylabel('phase-shift (s)','FontSize',18)
fmax = input('What is fmax? ');

% estimates of dtstar from the data
ind0 = min(find(fmids<=fmax));
fo1 = fit(fmids(ind0:end),log(As(ind0:end)),'poly1');
dtstar_e1 = -fo1.p1./pi;
fo2 = fit(log(fmids(ind0:end)),phis(ind0:end),'poly1');
dtstar_e2 = -fo2.p1*pi;
% estimates of dT from the data
dT_anel = fo2.p2 + fo2.p1*pi; % using phase lags, taking into account anelasticity
dT_xcor = diff(xcortimes([qdatwc1,qdatwc2], dt, 50, 10,0)); % ignoring anelasticity

figure(3),subplot(211),title(sprintf('True t* = %.3f  Est t* = %.3f',dtstar,dtstar_e1),'FontSize',20)
figure(3),subplot(212),title(sprintf('True t* = %.3f  Est t* = %.3f',dtstar,dtstar_e2),'FontSize',20)


% estimate from joint inversion of both
[ dtstar_e3,dT_e3,A0_e3 ] = invert_Aphi_4_dtdtstar( As(ind0:end),phis(ind0:end),fmids(ind0:end));

%% estimate dtstar from spectral ratio
wlen = length(qdatwc1); % window length, accounting for taper+padding
ntap=2;
nft=2^nextpow2(wlen);
[spec1,frq]=pmtm(qdatwc1,ntap,nft,samprate);
[spec2,~]=pmtm(qdatwc2,ntap,nft,samprate);
spec1 = spec1(2:length(spec1)).^0.5;
spec2 = spec2(2:length(spec2)).^0.5;
figure(4)
dtstar_lnR =  diff(xspecratio( [spec1,spec2],frq,fmax,0.01,0,1 ))

%% summarise
fprintf('\nPositive dt* means 2 more attenuated than 1\n')
fprintf('dt* true = %.3f\n',dtstar)
fprintf('dt* est1 = %.3f from amplitude\n',dtstar_e1)
fprintf('dt* est2 = %.3f from phase\n',dtstar_e2)
fprintf('dt* est3 = %.3f from both\n',dtstar_e3)
fprintf('dt* est4 = %.3f from spectral ratio\n\n',dtstar_lnR)

fprintf('\nPositive dT means 2 slower than 1\n')
fprintf('dT true\t= %.3f\n',dT)
fprintf('dT est1\t= %.3f with anelasticity phase \n',dT_anel)
fprintf('dT est2\t= %.3f with anelasticity both\n',dT_e3)
fprintf('dT est3\t= %.3f only elasticity\n',dT_xcor)

%% Use combspectra function
parms.comb.Tmin = 1; % fnq/2
parms.comb.Tmax = 20;
parms.comb.Nwds = 20;
parms.comb.Tw_opt = 'scale';
parms.comb.npol = 4;

parms.wind.pretime = 50;
parms.wind.prex = 5;
parms.wind.postx = 15;
parms.wind.taperx = 0.1;

parms.qc.minacor = 0.6;
parms.qc.maxphi = 8;

parms.inv.amp2phiwt = 1;
parms.inv.fmin = .2;
parms.inv.fmax = .7;
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
[delta_tstar,delta_T,std_dtstar,pairwise,frq] = combspectra_nofuss([qdat1,qdat2],samprate,parms,0);
fprintf('BETTER? combspectra dtstar est\t= %.3f \n',diff(delta_tstar))
fprintf('BETTER? combspectra dT est\t= %.3f \n',diff(delta_T))
toc



%% START THE FANCY PLOTS
figure(33), clf, set(gcf,'pos',[600 600 600,800])

% ================ PLOT THE WAVEFORMS ================
subplot(5,1,1), hold on
plot(tt,dat./max(abs(dat)),'k','LineWidth',2)
plot(tt,qdat1./max(abs(dat)),'b','LineWidth',2)
plot(tt,qdat2./max(abs(dat)),'r','LineWidth',2)
% plot(-5*[1 1],[0 1],'b--',15*[1 1],[0 1],'b--','linewidth',1.5)
text(-4,0.86,'\textbf{Original}','interpreter','latex','fontsize',14,'horizontalalignment','left','verticalalignment','bottom','color','k')
text(-4,0.63,sprintf('\\textbf{Q=%.0f}',Q0_1),'interpreter','latex','fontsize',14,'horizontalalignment','left','verticalalignment','bottom','color','b')
text(-4,0.4,sprintf('\\textbf{Q=%.0f}',Q0_2),'interpreter','latex','fontsize',14,'horizontalalignment','left','verticalalignment','bottom','color','r')
text(7,0.15,'Time $\rightarrow$','interpreter','latex','fontsize',15,'horizontalalignment','left','verticalalignment','bottom')
set(gca,'xlim',[-5 10],'color','none','ylim',[-0.1 1.1],...%'XAxisLocation','top',...
    'Fontsize',15,'linewidth',2,'box','on','ytick',[])

% ================ PLOT THE AMPLITUDE SPECTRA ================
subplot(5,1,2:3), hold on
plot(ff,log(A0_e3) - pi*ff*dtstar,'g','Linewidth',2.5)

scatter(frq,log(pairwise.As),150*pairwise.wts,'ok','MarkerFaceColor','r','linewidth',1.5)


% xlabel('freq (Hz)','FontSize',16,'interpreter','latex')
ylabel('$\ln\,(R_{12})$','FontSize',22,'interpreter','latex')
set(gca,'fontsize',15,'Xscale','linear','xlim',[0.00 1],'ylim',[-3 1],'linewidth',2,'box','on')

% ================ PLOT THE PHASE SPECTRA ================
subplot(5,1,4:5), hold on
scatter(frq,pairwise.phis,150*pairwise.wts,'ok','MarkerFaceColor','r','linewidth',1.5)
plot(ff,dtstar*(1 - log(ff)./pi) + dT,'g','Linewidth',2.5)

% ha1 = arrow([eqar(ind1).fcross,1],[eqar(ind1).fcross,0],5,90,2,0.1,'FaceColor','m'); % fcross for sta 1
% ha2 = arrow([eqar(ind2).fcross,1],[eqar(ind2).fcross,0],5,90,2,0.1,'FaceColor','g'); % fcross for sta 2
% plot(fmax*[1 1],[-1 1],'--b')  
text(0.025,-1.6,['$\mathbf{\Delta t^*_{obs} = ',num2str(diff(delta_tstar),2),'}$'],'interpreter','latex','fontsize',20,'horizontalalignment','left','verticalalignment','bottom')
text(0.025,-1,['$\mathbf{\Delta t^*_{tru} = ',num2str(dtstar,2),'}$'],'interpreter','latex','fontsize',20,'horizontalalignment','left','verticalalignment','bottom')

xlabel('frequency (Hz)','FontSize',22,'interpreter','latex')
ylabel('$\Delta \phi_{12}$ (s)','FontSize',22,'interpreter','latex')
set(gca,'fontsize',15,'Xscale','linear','xlim',[0.00 0.6],'ylim',[-2 3.5],'linewidth',2,'box','on')
return
% ================ SAVE FIGURE ================
ofile = sprintf('synth_comb_Q1-%.0f_Q2-%.0f',Q0_1,Q0_2);
save2pdf(33,ofile,'figs');


