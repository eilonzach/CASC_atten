% Second simple script to test the method of measuring differential
% attenuation using combs of filters and measuring amplitude loss as well
% as phase lag within each narrow bandpass.
% This script is attempts to do this with two traces that have passed
% through different Q and V material, to see how the velocity and Q
% differences trade off, and how the delta-tstar works.
% clear all
close all
cd /Users/zeilon/Documents/MATLAB/CASC_atten
addpath('matguts')
% parms
samprate = 40;
%% trace 1 Q & V - faster c0 and less attenuated
Q_1 = 200;
c0_1 = 3.9e3; % reference velocity in km/s
%% trace 2 Q & V - slower c0 and more attenuated
Q_2 = 20;
c0_2 = 3.6e3; % reference velocity in km/s

L = 100e3;

% calc. true values
dT = L*(1./c0_2 - 1./c0_1);
dtstar = L*(1./c0_2./Q_2 - 1./c0_1./Q_1);
fprintf('Theoretical dt=%.2f and dtstar=%.2f\n',dT,dtstar)
return
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
[ Dwt_1 ] = attenuation_operator( Q_1,c0_1,L,w);
[ Dwt_2 ] = attenuation_operator( Q_2,c0_2,L,w);
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
Nwds = 20;
Tmids = logspace(log10(2/fnq),log10(T/2),Nwds)';
Twdhs = 0.5*diff(logspace(log10(1/fnq),log10(T),Nwds+1)');
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
figure(5), clf, set(gcf,'pos',[600 10 500,700])
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

%% Test combspectra function
parms.comb.Tmin = 1; % fnq/2
parms.comb.Tmax = T/2;
parms.comb.Nwds = Nwds;
parms.comb.Tw_opt = 1;
parms.comb.npol = 4;

parms.wind.pretime = 50;
parms.wind.prex = 50;
parms.wind.postx = 50;
parms.wind.taperx = 0.1;

parms.qc.minacor = 0.5;
parms.qc.maxphi = 5;

parms.inv.amp2phiwt = 50;
parms.inv.fmin = 0.4;
parms.inv.fmax = 0.6;
parms.inv.ifwt = false;
parms.inv.corr_c_skip = true;

ifplot = false;

[delta_tstar,delta_T,~,pairwise] = combspectra_nofuss([qdat1,qdat2],samprate,parms,0);
fprintf('OLD combspectra dtstar est\t= %.3f \n',diff(delta_tstar))
fprintf('OLD combspectra dT est\t= %.3f \n',diff(delta_T))

% [delta_tstar_skip,delta_T_skip,~,pairwise_skip] = combspectra_cskip([qdat1,qdat2],samprate,parms,0);
% fprintf('BETTER? combspectra dtstar est\t= %.3f \n',diff(delta_tstar_skip))
% fprintf('BETTER? combspectra dT est\t= %.3f \n',diff(delta_T_skip))
% [delta_tstar_o,delta_T_o] = combspectra_old([qdat1,qdat2],samprate,parms,1);
