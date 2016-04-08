% Second simple script to test the method of measuring differential
% attenuation using combs of filters and measuring amplitude loss as well
% as phase lag within each narrow bandpass.
% This script is attempts to do this with five traces that have passed
% through different Q and V material, to see how the velocity and Q
% differences trade off, and how the delta-tstar works.
% clear all
% close all
addpath('matguts')
% parms
samprate = 20;
T = 2e3;

%% trace 1 Q & V - faster c0 and less attenuated
Q0_1 = 1000;
c0_1 = 3.8e3; % reference velocity in km/s
%% trace 2 Q & V - slower c0 and more attenuated
Q0_2 = 150;
c0_2 = 3.95e3; % reference velocity in km/s
%% trace 3 Q & V - faster c0 and less attenuated
Q0_3 = 50;
c0_3 = 3.9e3; % reference velocity in km/s
%% trace 4 Q & V - slower c0 and more attenuated
Q0_4 = 30;
c0_4 = 3.80e3; % reference velocity in km/s
%% trace 4 Q & V - slower c0 and more attenuated
Q0_5 = 20;
c0_5 = 3.80e3; % reference velocity in km/s

L = 200e3;

alpha = 0.3;


% calc. true values
TT = L*[c0_1,c0_2,c0_3,c0_4,c0_5].^-1;
dT_tru = demean(TT);

tstar = L*[c0_1*Q0_1,c0_2*Q0_2,c0_3*Q0_3,c0_4*Q0_4,c0_5*Q0_5].^-1;
dtstar_tru = demean(tstar);


%% ========== Make the two traces ==========

%% Make simple pulse
dt = 1/samprate;
tt = [-T/2:dt:T/2]'; % well padded
T = tt(end)-tt(1);
N = T/dt; tt = tt(1:N);
fnq = 0.5/dt;
dat = synthtrace(T,4,1,dt,'gauss');


%% Attenuate pulses
% take fft of pulse
[DAT,ff] = fft_ze(dat',dt);

w = 2*pi*ff;

% Make attenuation operators
[ Dwt_1 ] = attenuation_operator( Q0_1,c0_1,L,w,alpha,'zph');
[ Dwt_2 ] = attenuation_operator( Q0_2,c0_2,L,w,alpha,'zph');
[ Dwt_3 ] = attenuation_operator( Q0_3,c0_3,L,w,alpha,'zph');
[ Dwt_4 ] = attenuation_operator( Q0_4,c0_4,L,w,alpha,'zph');
[ Dwt_5 ] = attenuation_operator( Q0_5,c0_5,L,w,alpha,'zph');

% apply attenuation operator and ifft
qdat1 = real(ifft(DAT.*Dwt_1));
qdat2 = real(ifft(DAT.*Dwt_2));
qdat3 = real(ifft(DAT.*Dwt_3));
qdat4 = real(ifft(DAT.*Dwt_4));
qdat5 = real(ifft(DAT.*Dwt_5));

%% Delay pulses
qdat1 = interp1(tt+dT_tru(1),qdat1,tt,'linear',0); qdat1 = qdat1(:);
qdat2 = interp1(tt+dT_tru(2),qdat2,tt,'linear',0); qdat2 = qdat2(:);
qdat3 = interp1(tt+dT_tru(3),qdat3,tt,'linear',0); qdat3 = qdat3(:);
qdat4 = interp1(tt+dT_tru(4),qdat4,tt,'linear',0); qdat4 = qdat4(:);
qdat5 = interp1(tt+dT_tru(5),qdat5,tt,'linear',0); qdat5 = qdat5(:);

qdat = [qdat1,qdat2,qdat3,qdat4,qdat5];

%% Plot original and attenuated pulses
figure(1), clf, hold on
plot(tt,dat,'k','LineWidth',1.5)
plot(tt,qdat,'LineWidth',1.5)
title('Original (black), trace1 (blue) and trace2 (red)')
xlim([-10 10])

% %% Noisify
% addpath('~/Documents/MATLAB/seizmo-master/noise/')
% nspec = 10.^(4 + nlnm(abs(ff(2:end)))/20);
% SNR = 1000;
% qdat1 = real(ifft(DAT_1.*Dwt_1 .* SNR.*[1;nspec]));
% qdat2 = real(ifft(DAT_2.*Dwt_2 .* SNR.*[1;nspec]));
% figure(3), clf, hold on
% plot(tt,dat,'k','LineWidth',1.5)
% plot(tt,qdat1,'b','LineWidth',1.5)
% plot(tt,qdat2,'r','LineWidth',1.5)
% title('Original (black), trace1 (blue) and trace2 (red)')
% xlim([-10 10])


%% ========== do minimal example of combing ==========
parms.comb.Tmin = 1;
parms.comb.Tmax = 30;
parms.comb.Nwds = 30;
parms.comb.Tw_opt = 'scale';
parms.comb.npol = 2;

parms.wind.pretime = 1000;
parms.wind.prex = 10;
parms.wind.postx = 40;
parms.wind.taperx = 0.1;

parms.qc.minacor = 0.5;
parms.qc.maxphi = 15;

parms.inv.amp2phiwt = 1e3;
parms.inv.fmin = 0.05;
parms.inv.fmax = 1;
parms.inv.ifwt = true;
parms.inv.corr_c_skip = true;

parms.inv.alpha = alpha;

ifplot = true;

%% run the combs

[delta_tstar_1,delta_T_1,std_dtstar,pairwise,fmids] = combspectra_nofuss(qdat,samprate,parms,0);

Amat = pairwise.As;
phimat = pairwise.phis;
wtmat = double(pairwise.inds).*pairwise.wts;
test_alphas = [0:0.05:0.9];
% test_alphas = [0];
 
[ delta_tstar_2,delta_T_2,alpha_pref,alpha_misfits,dtstars,dTs,A0s, Eamp, Ephi ] ...
    = invert_Aphis_4_STA_dtdtstar_alpha(Amat,phimat,fmids,test_alphas,wtmat,parms.inv.amp2phiwt );

[dtstar_tru',delta_tstar_1,delta_tstar_2]
[dT_tru',delta_T_1,delta_T_2]
dtstars
dTs


