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
samprate = 20;
T = 2e3;
%% trace 1 Q & V - faster c0 and less attenuated
Q0_1 = 101;
c0_1 = 3.95e3; % reference velocity in km/s
%% trace 2 Q & V - slower c0 and more attenuated
Q0_2 = 20;
c0_2 = 3.80e3; % reference velocity in km/s

L = 200e3;

alpha = 0.27;

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

%% Plot original and attenuated pulses
figure(1), clf, hold on
plot(tt,dat,'k','LineWidth',1.5)
plot(tt,qdat1,'b','LineWidth',1.5)
plot(tt,qdat2,'r','LineWidth',1.5)
title('Original (black), trace1 (blue) and trace2 (red)')
xlim([-10 10])

%% Noisify
addpath('~/Documents/MATLAB/seizmo-master/noise/')
nspec = 10.^(4 + nlnm(abs(ff(2:end)))/20);
SNR = 1000;
qdat1 = real(ifft(DAT_1.*Dwt_1 .* SNR.*[1;nspec]));
qdat2 = real(ifft(DAT_2.*Dwt_2 .* SNR.*[1;nspec]));
figure(3), clf, hold on
plot(tt,dat,'k','LineWidth',1.5)
plot(tt,qdat1,'b','LineWidth',1.5)
plot(tt,qdat2,'r','LineWidth',1.5)
title('Original (black), trace1 (blue) and trace2 (red)')
xlim([-10 10])


%% ========== do minimal example of combing ==========
Tmin = 1;
Tmax = 30;
Nwds = 30;
npol = 2;

pretime = 1000;
prex    = 5;
postx   = 30;
taperx   = 0.1;

minacor = 0.7;
maxphi  = 10;

fmin      = 0.05;
fmax      = 0.5;
corc_skip = true;
ifwt      = true;

fnq = samprate/2;
dt = 1./samprate;
npt = size(qdat1,1);

ifplot = true;


%% prepare filter + cleaning parms
% Make set of period windows for bandpass filter
Tmids = logspace(log10(Tmin),log10(Tmax),Nwds)';
Twdhs = 0.5*diff(logspace(log10(Tmin/2),log10(2*Tmax),Nwds+1)');


% Twdhs = 0.5*Tmids;
fmids = 1./Tmids;

flos = 1./(Tmids + Twdhs);
fhis = 1./(Tmids - Twdhs);

clear fltinfo
for iw = 1:Nwds
    % option 1: butter
%     [fltdat(iw).bb,fltdat(iw).aa]=butter(parms.comb.npol, [flos(iw), fhis(iw)].*dt.*2);
    % option 2: cheby
%     [fltdat(iw).bb,fltdat(iw).aa]=cheby1(parms.comb.npol,0.5, [flos(iw), fhis(iw)].*dt.*2);
    % option 3: higher order butter
    [z,p,k]=butter(npol, [flos(iw), fhis(iw)].*dt.*2.);
    [sos,g]=zp2sos(z,p,k); 
    fltinfo(iw).bb=sos; 
    fltinfo(iw).aa=g;
    fltinfo(iw).fmid = fmids(iw);
    fltinfo(iw).flo = flos(iw);
    fltinfo(iw).fhi = fhis(iw);
end
% set up cleaning/taper/window parms
nwin=round((postx+prex)/dt); % window length in samples

wdo1 = tukeywin(npt,2*taperx);

n1=round((pretime-prex)/dt); % first sample in window
wdo2=[zeros(n1,1);tukeywin(nwin,2*taperx); zeros(npt-n1-nwin,1)]; % taperx% tukey window
ibds=[max(1,n1-floor(nwin*2*taperx)), min(npt,n1+nwin+floor(nwin*2*taperx))]; % extend so traces are extra 20% of window on either side
jbds=ibds(1):ibds(2); % indices of points to keep


[ As,phis,wts ] = run_comb( qdat1,qdat2,fltinfo,wdo1,wdo2,jbds,dt,pretime,maxphi,cor_c_skip,fmin,ifplot );

inds = find(fmids<=fmax & abs(phis)<maxphi & sqrt(wts)>minacor);

a_tests = [0:0.05:0.5]';
a_misfits = zeros(length(a_tests),1);
for ia = 1:length(a_tests)
[ ~,~,~,a_misfits(ia),~ ] = invert_Aphi_4_dtdtstar( As(inds),phis(inds),fmids(inds), wts(inds),3,a_tests(ia));
end
figure(7), clf, hold on
plot(a_tests,a_misfits)
fs = fit(a_tests,a_misfits,'smoothingspline');
plot(fs)
[ dtstar,dT] = invert_Aphi_4_dtdtstar( As(inds),phis(inds),fmids(inds), wts(inds),3,a_tests(mindex(a_misfits)))




%% =========== Use combspectra function ==================
parms.comb.Tmin = Tmin; % fnq/2
parms.comb.Tmax = Tmax;
parms.comb.Nwds = Nwds;
parms.comb.Tw_opt = 'scale';
parms.comb.npol = npol;

parms.wind.pretime = pretime;
parms.wind.prex = prex;
parms.wind.postx = postx;
parms.wind.taperx = taperx;

parms.qc.minacor = minacor;
parms.qc.maxphi = maxphi;

parms.inv.amp2phiwt = 3;
parms.inv.fmin = fmin;
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
[delta_tstar,delta_T,std_dtstar,pairwise,frq] = combspectra_nofuss([qdat1,qdat2],samprate,parms,0);
fprintf('Combspectra dtstar est\t= %.3f \n',diff(delta_tstar))
fprintf('Combspectra dT est\t= %.3f \n',diff(delta_T))
toc
return


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


