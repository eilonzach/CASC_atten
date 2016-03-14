% First simple script to test the method of measuring differential
% attenuation using combs of filters and measuring amplitude loss as well
% as phase lag within each narrow bandpass.
% This script is a proof of concept to measure tstar for one ideal,
% attenuated pulse
clear all
close all
% parms
samprate = 40;
Q = 100;
c0 = 5e3; % reference velocity in km/s
L = 300e3;

minf = 0.00001;
maxf = 5;
delf = 0.025;

%% ========== Simple pulse, check out phase delay ==========

%% Make pulse
dt = 1/samprate;
tt = [-50:dt:50]';
T = tt(end)-tt(1);
N = T/dt; tt = tt(1:N);
fnq = 0.5/dt;
dat = synthtrace(T,1,1,dt,'gauss');

%% Attenuate pulse
% take fft of pulse
[DAT,ff,DAT_0,ff_0] = fft_ze(dat,dt);
w = 2*pi*ff;
% Make attenuation operator
[ Dwt1 ] = attenuation_operator( Q,c0,L,w);
[ Dwt2 ] = attenuation_operator( 0.3*Q,0.98*c0,L,w);
% apply attenuation operator and ifft
qdat1 = ifft(DAT.*Dwt1);
qdat2 = ifft(DAT.*Dwt2);

%% calc. t-star
tstar = L./c0./Q;

%% Plot original and attenuated pulses
figure(1), clf, hold on
plot(tt,qdat1,'k','LineWidth',1.5)
plot(tt,qdat2,'b','LineWidth',1.5)


%% wavelet decomposition
inds = [1700:2499]';
Npts = length(inds);
s0 = 2*dt;
ds = 0.4875;
NbSc = fix(log2(Npts*dt/s0)/ds);
scales = s0*2.^((0:NbSc-1)*ds);
cwtSq1 = cwtft(qdat1(inds),'scales',scales,'wavelet','paul');
cwtSq2 = cwtft(qdat2(inds),'scales',scales,'wavelet','paul');

Nwds = length(cwtSq1.scales);

ttws = tt(inds);
re_cwtSq1 = cwtSq1; re_cwtSq1.cfs(:,:)=0;
re_cwtSq2 = cwtSq2; re_cwtSq2.cfs(:,:)=0;

fmids = cwtSq1.scales';
As = zeros(Nwds,1);
phis = zeros(Nwds,1);
for ii = 1:Nwds
    re_cwtSq1.cfs(:,:)=0;
    re_cwtSq2.cfs(:,:)=0;
    
    re_cwtSq1.cfs(1:ii,:) = cwtSq1.cfs(1:ii,:);
    re_cwtSq2.cfs(1:ii,:) = cwtSq2.cfs(1:ii,:);
    
    qdatwf1 = icwtft(re_cwtSq1)'; qdatwf1 = detrend(qdatwf1,'constant').*tukeywin(Npts,2*.1);
    qdatwf2 = icwtft(re_cwtSq2)'; qdatwf2 = detrend(qdatwf2,'constant').*tukeywin(Npts,2*.1);
        
    % find observed phase shift
    [dcor, dcstd, dvcstd, acor]=xcortimes([qdatwf1,qdatwf2], dt, 50, 10,0);
    phi_f_obs = diff(dcor);
    
    % make phase-corrected time series
    qdatwfs2 = interp1(ttws-phi_f_obs,qdatwf2,ttws,'linear',0);
    
    % calc. observed amplitude factor
    A_f_obs = (qdatwfs2'*qdatwf1)/(qdatwf1'*qdatwf1);
%     A_f_obs = 1./((datwf'*qdatwfs)/(datwf'*datwf)); %<< could do it other way too
    
    % make amplitude-corrected time series
    qdatwfsa2 = qdatwfs2./A_f_obs;

	figure(2),clf, hold on
    plot(ttws,qdatwf1,'k','LineWidth',2)
    plot(ttws,qdatwf2,'b','LineWidth',1.5)
    plot(ttws,qdatwfs2,'r','LineWidth',1.5)
    plot(ttws,qdatwfsa2,'g','LineWidth',1)
    pause

    As(ii) = A_f_obs;
    phis(ii) = phi_f_obs;
end

% figure(3), clf, set(gcf,'pos',[600 10 500,700])
% subplot(211), hold on
% plot(log(As),'or')
% xlabel('log(freq)','FontSize',18), ylabel('log(Amp)','FontSize',18)
% subplot(212), hold on
% plot(phis,'or')
% xlabel('log(freq)','FontSize',18), ylabel('phase-shift (s)','FontSize',18)
% 
% 
% 
% return
% %% ========== Load the two traces ==========
% % load('synthetic_tstar/eg_eqar_456ST.mat')
% load('synthetic_tstar/eg_eqar_456ST.mat')
% 
% ind1 = 6; %  1  1  1  2
% ind2 = 10;  % 48 -3  5 30
% 
% dat1 = eqar(ind1).datT';
% dat2 = eqar(ind2).datT';
% 
% inds = [7500:8500]';
% 
% cwtSd1= cwtft(dat1(inds),'scales',scales,'wavelet',wavelet,'plot');
% cwtSd2= cwtft(dat2(inds),'scales',scales,'wavelet',wavelet,'plot');
% figure(23), hold on
% plot(cwtSd1.cfs','LineWidth',2)
% plot(cwtSd2.cfs',':')
% 
% 
% %% Plot two traces
% figure(1), clf, set(gcf,'pos',[100 550 1000 350]), hold on
% plot(dat1./max(abs(dat1)),'k','LineWidth',1.5)
% plot(dat2./max(abs(dat2)),'r','LineWidth',1.5)
% 
% 
% return
% %% Make set of period windows for bandpass filter
% Nwds = 20;
% Tmids = logspace(log10(2/fnq),log10(T/2),Nwds)';
% Twdhs = 0.5*diff(logspace(log10(1/fnq),log10(T),Nwds+1)');
% fmids = 1./Tmids;
% As = zeros(Nwds,1);
% phis = zeros(Nwds,1);
% 
% for ii = 1:Nwds
%     flo = 1./(Tmids(ii) + Twdhs(ii));
%     fhi = 1./(Tmids(ii) - Twdhs(ii));
%     fmid = 1./Tmids(ii);
%     cp = struct('samprate',samprate,'pretime',50,'prex',50,'postx',50,...
%                 'taperx',0.1,'fhi',fhi,'flo',flo,'npoles',2,'norm',0);
%     [ datwf, datf, datwc, datc, ~, ttws, tts ] = data_clean( dat,cp );
%     [ qdatwf, qdatf, qdatwc, adatc, ~, ~, ~ ] = data_clean( qdat,cp );
%         
%     % find observed phase shift
%     [dcor, dcstd, dvcstd, acor]=xcortimes([datwf,qdatwf], dt, 50, 10,0);
%     phi_f_obs = diff(dcor);
%     % calc. predicted phase shift
%     phi_f_pred = tstar*log(fnq/fmid)./pi;
% 
%     
%     % make phase-corrected time series
%     qdatwfs = interp1(ttws-phi_f_obs,qdatwf,ttws,'linear',0)';
%     
%     % calc. observed amplitude factor
%     A_f_obs = (qdatwfs'*datwf)/(datwf'*datwf);
% %     A_f_obs = 1./((datwf'*qdatwfs)/(datwf'*datwf)); %<< could do it other way too
%     % calc. predicted amplitude factor from tstar
%     A_f_pred = exp(-pi*fmid*tstar);
%     
%     % make amplitude-corrected time series
%     qdatwfsa = qdatwfs./A_f_obs;
%     
%     figure(2), clf, set(gcf,'pos',[30 350 1200,400]), hold on
%     plot(ttws,datwf,'k','LineWidth',2)
%     plot(ttws,qdatwf,'b','LineWidth',1.5)
%     plot(ttws,qdatwfs,'r','LineWidth',1.5)
%     plot(ttws,qdatwfsa,'g','LineWidth',1)
%     
%     title(sprintf('Period = %.2f, A$_{pred}$ = %.2f, A$_{obs}$ = %.2f, $\\phi_{pred}$ = %.2f, $\\phi_{obs}$ = %.2f',...
%         Tmids(ii),A_f_pred,A_f_obs,phi_f_pred,phi_f_obs),'Fontsize',18,'interpreter','latex')
%     xlim([-20 20])
%     
%     As(ii) = A_f_obs;
%     phis(ii) = phi_f_obs;
% end

%% use vectors of amplitude factors and phase shifts to estimate tstar
% in theory, A(f)   = exp(pi*f*tstar)
%  so a plot of ln(A) against freq should have slope pi*tstar
%  or a plot of ln(A) against 1/t should have slope pi*tstar
% in theory, phi(f) = 1./pi * tstar * ln(fnq/f)
%  so a plot of phi against ln(freq) should have slope -tstar./pi
%  or a plot of phi against ln(t) should have slope tstar./pi
figure(3), clf, set(gcf,'pos',[600 10 500,700])
subplot(211), hold on
plot(fmids,log(As),'or')
plot(fmids,-pi*fmids*tstar,'b'), set(gca,'Xscale','log')
xlabel('log(freq)','FontSize',18), ylabel('log(Amp)','FontSize',18)
subplot(212), hold on
plot(fmids,phis,'or')
plot(fmids,log(fnq)*tstar./pi - log(fmids)*tstar./pi,'b'), set(gca,'Xscale','log')
xlabel('log(freq)','FontSize',18), ylabel('phase-shift (s)','FontSize',18)

fmax = input('What is fmax? ');

inds = fmids<=fmax;
fo1 = fit(fmids(inds),log(As(inds)),'poly1');
tstar_e1 = -fo1.p1./pi;
fo2 = fit(log(fmids(inds)),phis(inds),'poly1');
tstar_e2 = -fo2.p1*pi;
fprintf('True t* = %.3f\nEst1 t* = %.3f\nEst2 t* = %.3f\n',tstar,tstar_e1,tstar_e2)
figure(3),subplot(211),title(sprintf('True t* = %.3f  Est t* = %.3f',tstar,tstar_e1),'FontSize',20)
figure(3),subplot(212),title(sprintf('True t* = %.3f  Est t* = %.3f',tstar,tstar_e2),'FontSize',20)
