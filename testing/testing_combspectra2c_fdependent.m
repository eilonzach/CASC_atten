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
samprate = 50;
T = 2e3;
%% trace 1 Q & V - faster c0 and less attenuated
Q0_1 = 101;
c0_1 = 3.95e3; % reference velocity in km/s
%% trace 2 Q & V - slower c0 and more attenuated
Q0_2 = 20;
c0_2 = 3.80e3; % reference velocity in km/s

L = 200e3;

alpha = 0.6;

fmax = .2;

xx = 0.5*cot(alpha*pi/2)*(2*pi).^alpha;

% calc. true values
dT = L*(1./c0_2 - 1./c0_1);
dtstar = L*(1./c0_2./Q0_2 - 1./c0_1./Q0_1);
fprintf('Theoretical dT0=%.2f and dtstar0=%.2f\n',dT,dtstar)

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
[DAT_1,~ ] = fft_ze(dat',dt);
[DAT_2,ff] = fft_ze(dat',dt);
w = 2*pi*ff;
% Make attenuation operators
[ Dwt_1 ] = attenuation_operator( Q0_1,c0_1,L,w,alpha,'zph');
[ Dwt_2 ] = attenuation_operator( Q0_2,c0_2,L,w,alpha,'zph');
% apply attenuation operator and ifft
qdat1 = real(ifft(DAT_1.*Dwt_1));
qdat2 = real(ifft(DAT_2.*Dwt_2));

%% Delay pulses
dT1 = L*(1./c0_1 - 1./mean([c0_1,c0_2]));
dT2 = L*(1./c0_2 - 1./mean([c0_1,c0_2]));
qdat1 = interp1(tt,qdat1,tt,'linear',0);           qdat1 = qdat1(:);
qdat2 = interp1(tt+(dT2-dT1),qdat2,tt,'linear',0); qdat2= qdat2(:);


% %% plot and test attenuation operators' relative phase and amplitude
% figure(21), clf, set(gcf,'position',[-20    21   999   492])
% 
% subplot(211), hold on
% plot(w(1:N/2).^(1-alpha),log(abs(Dwt_1(1:N/2))),'b')
% plot(w(1:N/2).^(1-alpha),log(abs(Dwt_2(1:N/2))),'r')
% RR = Dwt_2./Dwt_1;
% plot(w(1:N/2).^(1-alpha),log(abs(RR(1:N/2))),'g')
% plot(w(1:N/2).^(1-alpha),-.5*w(1:N/2).^(1-alpha)*dtstar,'go')
% set(gca,'xscale','linear','xlim',[.0 1])
% 
% subplot(212),  hold on
% semilogx(w(1:N/2).^(-alpha),angle(Dwt_1(1:N/2)),'b')
% semilogx(w(1:N/2).^(-alpha),angle(Dwt_2(1:N/2)),'r')
% plot(w(1:N/2).^(-alpha),angle(Dwt_2(1:N/2)./Dwt_1(1:N/2)),'g')
% plot(w(1:N/2).^(-alpha),mp2pp(-dtstar*xx*w(1:N/2).^(1-alpha)),'go')
% set(gca,'xscale','log','xlim',[.05 100])
% 


%% Plot original and attenuated pulses
figure(1), clf, hold on
plot(tt,dat,'k','LineWidth',1.5)
plot(tt,qdat1,'b','LineWidth',1.5)
plot(tt,qdat2,'r','LineWidth',1.5)
title('Original (black), trace1 (blue) and trace2 (red)')
xlim([-10 10])

%% Make set of period windows for bandpass filter
Tmax = T/40;
Tmin = 1;
Nwds = 20;
Tmids = logspace(log10(Tmin),log10(Tmax),Nwds)';
Twdhs = 0.5*diff(logspace(log10(Tmin),log10(Tmax),Nwds+1)');
fmids = 1./Tmids;
wmids = 2*pi*fmids;
As = zeros(Nwds,1);
phis = zeros(Nwds,1);

for ii = 1:Nwds
    flo = 1./(Tmids(ii) + Twdhs(ii));
    fhi = 1./(Tmids(ii) - Twdhs(ii));
    fmid = 1./Tmids(ii);
    wmid = 2*pi*fmid;
    cp = struct('samprate',samprate,'pretime',T/2,'prex',T/2,'postx',T/2,...
                'taperx',0.1,'fhi',fhi,'flo',flo,'npoles',2,'norm',0);
    [ qdatwf1, qdatf1, qdatwc1, ~, ~, ttws, ~ ] = data_clean( qdat1,cp );
    [ qdatwf2, qdatf2, qdatwc2, ~, ~, ~,    ~ ] = data_clean( qdat2,cp );

    % find observed phase shift
    [dcor, dcstd, dvcstd, acor]=xcortimes([qdatwf1,qdatwf2], dt, T/2, 10,0);
    phi_f_obs = diff(dcor);
    % calc. predicted phase shift
%     phi_f_pred = dtstar*log(fnq/fmid)./pi + dT; % alpha==0
    phi_f_pred = dtstar*xx*wmid.^(1-alpha) + dT;
    phi_f_pred = dtstar*xx*wmid.^(-alpha) + dT;

    
    % make phase-corrected time series
    qdatwf2s = interp1(ttws-phi_f_obs,qdatwf2,ttws,'linear',0)';
    % calc. observed amplitude factor
    A_f_obs = (qdatwf2s'*qdatwf1)/(qdatwf1'*qdatwf1);
%     A_f_obs = 1./((qdatwf1'*qdatwf2s)/(qdatwf2s'*qdatwf2s)); %<< could do it other way too
    % calc. predicted amplitude factor from tstar
    A_f_pred = exp(-.5*wmid.^(1-alpha)*dtstar);

    % make amplitude-corrected time series
    qdatwf2sa = qdatwf2s./A_f_obs;
    

    As(ii) = A_f_obs;
    phis(ii) = phi_f_obs;
    
%     if mod(ii,2) % only plot half the cfs
%     figure(2), clf, set(gcf,'pos',[30 350 1200,400]), hold on
%     plot(ttws,qdatwf1,'b','LineWidth',2)
%     plot(ttws,qdatwf2,'r','LineWidth',1.5)
%     plot(ttws,qdatwf2s,'m','LineWidth',1.5)
%     plot(ttws,qdatwf2sa,'g','LineWidth',1)
%     legend('dat1','dat2','dat2s','dat2sa')
%     
%     title(sprintf('Period = %.2f, A$_{pred}$ = %.2f, A$_{obs}$ = %.2f, $\\phi_{pred}$ = %.2f, $\\phi_{obs}$ = %.2f',...
%         Tmids(ii),A_f_pred,A_f_obs,phi_f_pred,phi_f_obs),'Fontsize',18,'interpreter','latex')
%     xlim([-100 100])
% %      pause
%     end
end


%% use vectors of amplitude factors and phase shifts to estimate dtstar and dT
% in theory, A(f)   = exp(pi*f*tstar)
%  so a plot of ln(A) against freq should have slope pi*dtstar
%  or a plot of ln(A) against 1/t should have slope pi*dtstar
% in theory, phi(f) = 1./pi * tstar * ln(fnq/f)
%  so a plot of phi against ln(freq) should have slope -dtstar./pi and
%  intercept that related to fnq and dt
%  or a plot of phi against ln(t) should have slope dtstar./pi

ff = logspace(log10(fmids(end)),log10(fmids(1)),40);
ww = 2*pi*ff;

%% plot spectra and fits
figure(3), clf, set(gcf,'pos',[600 10 500,700])
subplot(211), hold on
plot(fmids.^(1-alpha),log(As),'or')
plot(ff.^(1-alpha),-0.5*(2*pi*ff).^(1-alpha).*dtstar,'r--'), 

set(gca,'Xscale','linear')
title(num2str(alpha))
xlabel('freq\^(1-a)','FontSize',18), ylabel('log(Amp)','FontSize',18)

subplot(212), hold on
plot(fmids.^(-alpha),phis,'or')
plot(ff.^(-alpha),dtstar*xx*(2*pi*ff).^(-alpha) + dT,'r--'), 
set(gca,'Xscale','linear')
xlabel('freq\^(-a)','FontSize',18), ylabel('phase-shift (s)','FontSize',18)


% subplot(212), hold on
% plot(ff,1.18*dtstar.*xx.*(2*pi).^(-alpha)*ff.^(-0.8*alpha),'g--'), 
% plot(ff,((2*pi).^(-alpha)*ff.^(-alpha)).*(L*(1./c0_2/Q0_2 - 1./c0_1/Q0_1).*xx),'g--'), 
% plot(ff.^(-alpha), (2*pi).^(-alpha)*ff.^(-alpha).*(...
%      xx*dtstar...
%      ),'c--'), 
plot(ff.^(-alpha), (2*pi).^(-alpha)*ff.^(-alpha).*(...
     xx*dtstar...
     + xx*xx*L/c0_2/Q0_2/Q0_2 - xx*xx*L/c0_1/Q0_1/Q0_1 ...
     ) + dT,'m--'), 
plot(ff.^(-alpha), (2*pi).^(-alpha)*ff.^(-alpha).*(...
     xx*dtstar...
     + xx*xx*L/c0_2/Q0_2/Q0_2 - xx*xx*L/c0_1/Q0_1/Q0_1 ...
     + xx*xx*xx*L/c0_2/Q0_2/Q0_2/Q0_2 - xx*xx*xx*L/c0_1/Q0_1/Q0_1/Q0_1...
     ) + dT,'g--'), 
% estimates of dtstar from the data
ind0 = min(find(fmids<=fmax)); %#ok<MXFND>
fo1 = fit(fmids(ind0:end).^(1-alpha),log(As(ind0:end)),'poly1');
dtstar_e1 = -fo1.p1./(2.^(-alpha)*pi.^(1-alpha))

%% test values of alpha
alpha_tests = [0.05:0.05:0.5]';
dtstar_est = zeros(length(alpha_tests),1);
misfit = zeros(length(alpha_tests),1);
for ia = 1:length(alpha_tests)
[foa,gda] = fit(fmids(ind0:end).^(1-alpha_tests(ia)),log(As(ind0:end)),'poly1');
dtstar_est(ia) = foa.p1./(-0.5*(2*pi).^(1-alpha_tests(ia)));
figure(4),clf, hold on
plot(fmids, log(As), 'ok')
plot(fmids(ind0:end), log(As(ind0:end)), 'ok', 'MarkerFaceColor','r')
plot(fmids,foa.p1*fmids.^(1-alpha) + foa.p2,'k--')
% plot(fmids,-0.5*(2*pi*ff).^(1-alpha_tests(ia)).*dtstar
misfit(ia) = gda.sse;
end    
alpha_e1 = alpha_tests(mindex(misfit));
dtstarae_e1 = dtstar_est(mindex(misfit))

figure(5)
plot(alpha_tests,misfit)



fo2 = fit(fmids(ind0:end).^(-alpha),phis(ind0:end),'poly1');
dtstar_e2 = fo2.p1./xx./(2*pi).^(-alpha)
% estimates of dT from the data
dT_anel = fo2.p2/(2*pi).^(-alpha); % using phase lags, taking into account anelasticity
dT_xcor = diff(xcortimes([qdatwc1,qdatwc2], dt, T/2, 10,0)); % ignoring anelasticity



% estimate from joint inversion of both
% [ dtstar_e3,dT_e3,A0_e3 ] = invert_Aphi_4_dtdtstar( As(ind0:end),phis(ind0:end),fmids(ind0:end))
[ dtstar_e3,dT_e3,A0_e3,misfit, E ] =...
    invert_Aphi_4_dtdtstar_testing( As(ind0:end),phis(ind0:end),fmids(ind0:end),[],5,.25);
dtstar_e3


figure(3)
subplot(211), hold on
plot(ff.^(1-alpha),-2.^(-alpha)*pi.^(1-alpha)*dtstar_e1*ff.^(1-alpha),'b')
plot(ff.^(1-alpha),-2.^(-alpha)*pi.^(1-alpha)*dtstarae_e1*ff.^(1-alpha),'b--')
plot(ff.^(1-alpha),-2.^(-alpha)*pi.^(1-alpha)*dtstar_e2*ff.^(1-alpha),'k')
plot(ff.^(1-alpha),-2.^(-alpha)*pi.^(1-alpha)*dtstar_e3*ff.^(1-alpha) + log(A0_e3),'g')
set(gca,'Xscale','linear')

subplot(212), hold on
plot(ff.^(-alpha),2*ff.^(-alpha)*dtstar_e2*((2*pi).^(alpha-1))*cot(alpha*pi/2) + dT_anel,'k')
plot(ff.^(-alpha),2*ff.^(-alpha)*dtstar_e3*((2*pi).^(alpha-1))*cot(alpha*pi/2) + dT_e3,'g')
set(gca,'Xscale','linear')

figure(3),subplot(211),title(sprintf('True t* = %.3f  Est t* = %.3f',dtstar,dtstar_e3),'FontSize',20)
figure(3),subplot(212),title(sprintf('True dT = %.3f  Est dT = %.3f',dT,dT_e3),'FontSize',20)


%% estimate dtstar from spectral ratio
wlen = length(qdatwc1); % window length, accounting for taper+padding
ntap=2;
nft=2^nextpow2(wlen);
[spec1,frq]=pmtm(qdatwc1,ntap,nft,samprate);
[spec2,~]=pmtm(qdatwc2,ntap,nft,samprate);
spec1 = spec1(2:length(spec1)).^0.5;
spec2 = spec2(2:length(spec2)).^0.5;
figure(4)
dtstar_lnR =  diff(xspecratio( [spec1,spec2],frq,fmax,0.01,0,1 ))'

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


