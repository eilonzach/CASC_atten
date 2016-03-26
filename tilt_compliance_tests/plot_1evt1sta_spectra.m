function plot_1evt1sta_spectra(spectra_file)
% code run after OBS_noise_cohere that reviews and plots spectral
% properties for a single station
% 
% INPUTS
%   spectra_file = .mat file with all spectra details for one station for
%                   one event, e.g. % spectra_file = '/Volumes/DATA/CASCADIA/DATA/215_201308301625/J32C_spectra.mat';
close all

load(spectra_file);

grey = [0.5 0.5 0.5];


% Plotting Power
figure(3), clf, set(gcf,'pos',[30 30 800 600])
subplot(411);
loglog(f,smooth(czz_stack,40),'-','LineWidth',.5,'Color',grey)
title('Power','Fontsize',16,'Interpreter','latex')
subplot(412);
loglog(f,smooth(c11_stack,40),'-','LineWidth',.5,'Color',grey)
subplot(413);
loglog(f,smooth(c22_stack,40),'-','LineWidth',.5,'Color',grey)
subplot(414);
loglog(f,smooth(cpp_stack,40),'-','LineWidth',.5,'Color',grey)


% Plotting Coherence
figure(4), clf, set(gcf,'pos',[40 40 800 600])
subplot(611);
semilogx(f,smooth(coh1z_stack,40),'-','LineWidth',.5,'Color',grey);
title('Coherence','Fontsize',16,'Interpreter','latex')
subplot(612);
semilogx(f,smooth(coh2z_stack,40),'-','LineWidth',.5,'Color',grey);
subplot(613);
semilogx(f,smooth(cohpz_stack,40),'-','LineWidth',.5,'Color',grey);
subplot(614);
semilogx(f,smooth(coh12_stack,40),'-','LineWidth',.5,'Color',grey);
subplot(615);
semilogx(f,smooth(coh1p_stack,40),'-','LineWidth',.5,'Color',grey);
subplot(616);
semilogx(f,smooth(coh2p_stack,40),'-','LineWidth',.5,'Color',grey);


% Plotting Phase
figure(6), clf, set(gcf,'pos',[60 60 800 600])
subplot(611);
semilogx(f,ph1z_stack,'-','LineWidth',.5,'Color',grey);
title('Phase','Fontsize',16,'Interpreter','latex')
subplot(612);
semilogx(f,ph2z_stack,'-','LineWidth',.5,'Color',grey);
subplot(613);
semilogx(f,phpz_stack,'-','LineWidth',.5,'Color',grey);
subplot(614);
semilogx(f,ph12_stack,'-','LineWidth',.5,'Color',grey);
subplot(615);
semilogx(f,ph1p_stack,'-','LineWidth',.5,'Color',grey);
subplot(616);
semilogx(f,ph2p_stack,'-','LineWidth',.5,'Color',grey);


% Plotting Admittance
figure(5), clf, set(gcf,'pos',[50 50 800 600])
subplot(611);
loglog(f,smooth(ad1z_stack,40),'-','LineWidth',.5,'Color',grey);
title('Admittance','Fontsize',16,'Interpreter','latex')
subplot(612);
loglog(f,smooth(ad2z_stack,40),'-','LineWidth',.5,'Color',grey);
subplot(613);
loglog(f,smooth(adpz_stack,40),'-','LineWidth',.5,'Color',grey);
subplot(614);
loglog(f,smooth(ad12_stack,40),'-','LineWidth',.5,'Color',grey);
subplot(615);
loglog(f,smooth(ad1p_stack,40),'-','LineWidth',.5,'Color',grey);
subplot(616);
loglog(f,smooth(ad2p_stack,40),'-','LineWidth',.5,'Color',grey);

