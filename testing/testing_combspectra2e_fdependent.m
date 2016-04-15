% Second simple script to test the method of measuring differential
% attenuation using combs of filters and measuring amplitude loss as well
% as phase lag within each narrow bandpass.
% This script is attempts to do this with five traces that have passed
% through different Q and V material, to see how the velocity and Q
% differences trade off, and how the delta-tstar works.
clear all
close all
addpath('../matguts')
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
c0_4 = 3.85e3; % reference velocity in km/s
%% trace 4 Q & V - slower c0 and more attenuated
Q0_5 = 60;
c0_5 = 3.80e3; % reference velocity in km/s

L = 200e3;

alpha = 0.09;


% calc. true values
TT = L*[c0_1,c0_2,c0_3,c0_4,c0_5].^-1;
dT_tru = demean(TT)';

tstar = L*[c0_1*Q0_1,c0_2*Q0_2,c0_3*Q0_3,c0_4*Q0_4,c0_5*Q0_5].^-1;
dtstar_tru = demean(tstar)';


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

parms.inv.amp2phiwt = 1;
parms.inv.fmin = 0.05;
parms.inv.fmax = 1;
parms.inv.ifwt = true;
parms.inv.corr_c_skip = true;

parms.inv.alpha = alpha;

ifplot = true;

%% run the combs

[delta_tstar_1,delta_T_1,std_dtstar,pairwise,fmids] = combspectra(qdat,samprate,parms,0);

Amat = pairwise.As;
phimat = pairwise.phis;
wtmat = double(pairwise.inds).*pairwise.wts;
test_alphas = [0:0.05:0.4];
% test_alphas = [0.1];
 
%% all in one method
[ delta_tstar_1,delta_T_1,A0_1,alpha_pref_1,alpha_misfits,dtstars,dTs,A0s, Eamp, Ephi ] ...
    = invert_allin1_Aphis_4_STA_dtdtstar_alpha(Amat,phimat,fmids,test_alphas,wtmat,parms.inv.amp2phiwt );

% report the best fitting vs. estimated values
[dtstar_tru,delta_tstar_1,nan(size(A0_1)),dT_tru,delta_T_1, nan(size(A0_1)),A0_1]

for ia = 1:length(test_alphas)
% compute A and phi predictions for best fitting values at each alpha
[ Amat_pred1,phimat_pred1 ] = pred_Amat_phimat( dtstars(:,ia),dTs(:,ia),A0s(:,ia),fmids,test_alphas(ia) );
Eaw1(ia) = sum(sum((log(Amat_pred1)-log(Amat)).^2.*wtmat));
Epw1(ia) = sum(sum((phimat_pred1-phimat).^2.*wtmat));
end
figure(78)
plot(test_alphas,[Eaw1;Epw1])

%% one-by-one method
[ delta_tstar_2,delta_T_2,A0_2,alpha_pref_2,alpha_misfits,dtstars,dTs,A0s, Eamp, Ephi ] ...
    = invert_1by1_Aphis_4_STA_dtdtstar_alpha(Amat,phimat,fmids,test_alphas,wtmat,parms.inv.amp2phiwt );

% report the best fitting vs. estimated values
fprintf('dt*_tru  dt*_est  dT_tru  dT_est  A0_est   \n')
for ii = 1:length(dtstar_tru)
    fprintf(' %5.2f   %5.2f    %5.2f   %5.2f   %5.2f \n',...
            dtstar_tru(ii),delta_tstar_2(ii),dT_tru(ii),delta_T_2(ii),A0_2(ii))
end

for ia = 1:length(test_alphas)
% compute A and phi predictions for best fitting values at each alpha
[ Amat_pred2,phimat_pred2 ] = pred_Amat_phimat( dtstars(:,ia),dTs(:,ia),A0s(:,ia),fmids,test_alphas(ia) );
Eaw2(ia) = sum(sum((log(Amat_pred2)-log(Amat)).^2.*wtmat));
Epw2(ia) = sum(sum((phimat_pred2-phimat).^2.*wtmat));
end
figure(79)
plot(test_alphas,[Eaw2;Epw2])

%% predict matrices to assess misfit
[ Amat_pred,phimat_pred ] = pred_Amat_phimat( dtstar_tru,dT_tru,ones(size(dT_tru)),fmids,alpha );



%% Is misfit difference significant?!
ia = find(test_alphas==alpha_pref);

[ dtstar,dT,A0] = invert_1pair_Aphi_4_dtdtstar( Amat(4,:)',phimat(4,:)',fmids,[],1,alpha)

[~,~,Ea2,Ep2] = plot_AmpPhi_fit(Amat(4,:)',phimat(4,:)',fmids,[],...
    dtstar,...
    dT,...
    A0,...
    alpha,30);

[~,~,Ea1,Ep1] = plot_AmpPhi_fit(Amat(4,:)',phimat(4,:)',fmids,[],...
    dtstar_tru(5)-dtstar_tru(1),...
    dT_tru(5)-dT_tru(1),...
    1,...
    alpha,32);

res1 = [parms.inv.amp2phiwt*Ea1;Ep1].*[wtmat(1,:)';wtmat(1,:)'].^0.5;
res2 = [parms.inv.amp2phiwt*Ea2;Ep2].*[wtmat(1,:)';wtmat(1,:)'].^0.5;
ftest(res2,3,res1,3)
% yes.

%% compare all-in-1 to 1-by-1
[ dtstar_,dT_,A0_,alpha_,chi2,a,b,c,mf_a1,mf_p1] = invert_1by1_Aphis_4_STA_dtdtstar_alpha( Amat,phimat,fmids,[0:0.05:0.4],[],5)
[ dtstar_,dT_,A0_,alpha_,chi2,a,b,c,mf_a2,mf_p2] = invert_allin1_Aphis_4_STA_dtdtstar_alpha( Amat,phimat,fmids,[0:0.05:0.4],[],5)


[~,~,Ea2,Ep2] = plot_AmpPhi_fit(Amat(4,:)',phimat(4,:)',fmids,[],...
    dtstar_(5)-dtstar_(1),...
    dT_(5)-dT_(1),...
    A0_(5)/A0_(1),...
    alpha_,31);


for ia = 1:5
plot_AmpPhi_fit(Amat(2,:)',phimat(2,:)',fmids,[],...
    a(3,ia)-a(1,ia),...
    b(3,ia)-b(1,ia),...
    c(3,ia)/c(1,ia),...
    test_alphas(ia),33+ia);
end
return




