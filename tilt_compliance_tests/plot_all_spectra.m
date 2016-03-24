%  plot_spec
% code run after OBS_noise_cohere that reviews and plots spectral
% properties for a single station

clear all;
close all


infile = '/Volumes/DATA/CASCADIA/DATA/215_201308301625/J32C_spectra.mat'
figoutpath=sprintf('/Users/zeilon/Documents/MATLAB/CASC_atten/figs/spectra/%s_spectra',station);
outpath = sprintf('/Users/helenj/Cascadia/EARTHQUAKES/NOISETC/AVG_STA/');

if ~exist(figoutpath)
    mkdir(figoutpath);
end

if ~exist(outpath)
    mkdir(outpath);
end

cc = [0.5 0.5 0.5];
spectra_filenames = dir(fullfile(infile,['*.mat']));

for ie = 1 : length(spectra_filenames)
    % for ie = 1:1
    load(ifile);
    freqcomp = sqrt(9.8/(2*pi*elev));
    a = length(station);
    eventid = spectra_filenames(ie).name((a+2):(a+13));
        disp(eventid);

    
    maxpowz(ie) = max(czz_stack)*10;
    maxpow1(ie) = max(c11_stack)*10;
    maxpow2(ie) = max(c22_stack)*10;
    maxpowp(ie) = max(cpp_stack)*10;
    
    minpowz(ie) = min(czz_stack)/10;
    minpow1(ie) = min(c11_stack)/10;
    minpow2(ie) = min(c22_stack)/10;
    minpowp(ie) = min(cpp_stack)/10;
    
    % Plotting Power
    figure(3)
    subplot(411);
    loglog(f,smooth(czz_stack,40),'-','LineWidth',.5,'Color',cc)
    hold on
    subplot(412);
    loglog(f,smooth(c11_stack,40),'-','LineWidth',.5,'Color',cc)
    hold on
    subplot(413);
    loglog(f,smooth(c22_stack,40),'-','LineWidth',.5,'Color',cc)
    hold on
    subplot(414);
    loglog(f,smooth(cpp_stack,40),'-','LineWidth',.5,'Color',cc)
    hold on
    
    czz_all(ie,:) = czz_stack;
    cpp_all(ie,:) = cpp_stack;
    c11_all(ie,:) = c11_stack;
    c22_all(ie,:) = c22_stack;
    
    if check_spec ==1;
    in1 = input('Would you like to use data from this day (n/y)? ','s');
    if in1 == 'n';
        in2 = input('Are you sure you want to delete this day (y)? ','s');
        if in2 == 'y'
            delete(fullfile(infile,spectra_filenames(ie).name));
            continue
        end
    end
    end
    
    % Plotting Coherence
    figure(4)
    subplot(611);
    semilogx(f,smooth(coh1z_stack,40),'-','LineWidth',.5,'Color',cc);
    hold on
    subplot(612);
    semilogx(f,smooth(coh2z_stack,40),'-','LineWidth',.5,'Color',cc);
    hold on
    subplot(613);
    semilogx(f,smooth(cohpz_stack,40),'-','LineWidth',.5,'Color',cc);
    hold on
    subplot(614);
    semilogx(f,smooth(coh12_stack,40),'-','LineWidth',.5,'Color',cc);
    hold on
    subplot(615);
    semilogx(f,smooth(coh1p_stack,40),'-','LineWidth',.5,'Color',cc);
    hold on
    subplot(616);
    semilogx(f,smooth(coh2p_stack,40),'-','LineWidth',.5,'Color',cc);
    hold on
    
    coh1z_all(ie,:) = coh1z_stack;
    coh2z_all(ie,:) = coh2z_stack;
    cohpz_all(ie,:) = cohpz_stack;
    coh12_all(ie,:) = coh12_stack;
    coh1p_all(ie,:) = coh1p_stack;
    coh2p_all(ie,:) = coh2p_stack;
    
    if check_spec ==1
    in1 = input('Would you like to use data from this day ((a)ll/y)? ','s');
    if in1 == 'n';
        in2 = input('Are you sure you want to delete this day (y)? ','s');
        if in2 == 'y'
            delete(fullfile(infile,spectra_filenames(ie).name));
            continue
        end
    end
    end
    
    % Plotting Phase
    figure(6)
    subplot(611);
    semilogx(f,ph1z_stack,'-','LineWidth',.5,'Color',cc);
    hold on
    subplot(612);
    semilogx(f,ph2z_stack,'-','LineWidth',.5,'Color',cc);
    hold on
    subplot(613);
    semilogx(f,phpz_stack,'-','LineWidth',.5,'Color',cc);
    hold on
    subplot(614);
    semilogx(f,ph12_stack,'-','LineWidth',.5,'Color',cc);
    hold on
    subplot(615);
    semilogx(f,ph1p_stack,'-','LineWidth',.5,'Color',cc);
    hold on
    subplot(616);
    semilogx(f,ph2p_stack,'-','LineWidth',.5,'Color',cc);
    hold on
    
    
    % Plotting Admittance
    figure(5)
    subplot(611);
    loglog(f,smooth(ad1z_stack,40),'-','LineWidth',.5,'Color',cc);
    hold on
    subplot(612);
    loglog(f,smooth(ad2z_stack,40),'-','LineWidth',.5,'Color',cc);
    hold on
    subplot(613);
    loglog(f,smooth(adpz_stack,40),'-','LineWidth',.5,'Color',cc);
    hold on
    subplot(614);
    loglog(f,smooth(ad12_stack,40),'-','LineWidth',.5,'Color',cc);
    hold on
    subplot(615);
    loglog(f,smooth(ad1p_stack,40),'-','LineWidth',.5,'Color',cc);
    hold on
    subplot(616);
    loglog(f,smooth(ad2p_stack,40),'-','LineWidth',.5,'Color',cc);
    hold on
    
end

% Printing and label properties
figure(3)
subplot(411);
title(sprintf('Z-component, Station: %s',station));
xlabel('Frequency (Hz)')
ylabel('Power Spectrum')
ylim([10^-20 max(maxpowz)])
plot([pper pper],[min(minpowz) max(maxpowz)],'-b','LineWidth',2)
hold on
plot([sper sper],[min(minpowz) max(maxpowz)],'r','LineWidth',2)
hold on
plot([freqcomp freqcomp],[min(minpowz) max(maxpowz)],'-k','LineWidth',2)
subplot(412);
title(sprintf('H1-component, Station: %s',station));
xlabel('Frequency (Hz)')
ylabel('Power Spectrum')
ylim([10^-20 max(maxpow1)])
plot([pper pper],[min(minpow1) max(maxpowz)],'-b','LineWidth',2)
hold on
plot([sper sper],[min(minpow1) max(maxpowz)],'r','LineWidth',2)
hold on
plot([freqcomp freqcomp],[min(minpow1) max(maxpowz)],'-k','LineWidth',2)
subplot(413);
title(sprintf('H2-component, Station: %s',station));
xlabel('Frequency (Hz)')
ylabel('Power Spectrum')
ylim([10^-20 max(maxpow2)])
plot([pper pper],[min(minpow2) max(maxpowz)],'-b','LineWidth',2)
hold on
plot([sper sper],[min(minpow2) max(maxpowz)],'r','LineWidth',2)
hold on
plot([freqcomp freqcomp],[min(minpow2) max(maxpowz)],'-k','LineWidth',2)
subplot(414);
title(sprintf('P-component, Station: %s',station));
xlabel('Frequency (Hz)')
ylabel('Power Spectrum')
ylim([10^-10 max(maxpowp)])
plot([pper pper],[min(minpowp) max(maxpowp)],'-b','LineWidth',2)
hold on
plot([sper sper],[min(minpowp) max(maxpowp)],'r','LineWidth',2)
hold on
plot([freqcomp freqcomp],[min(minpowp) max(maxpowp)],'-k','LineWidth',2)
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);
filename=sprintf('%s/Spectra_%s',figoutpath,station);
print(gcf,'-dpng',filename)

figure(4)
subplot(611);
title(sprintf('Z-1, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Coherence')
ylim([0 1])
plot([pper pper],[0 1],'-b','LineWidth',2)
hold on
plot([sper sper],[0 1],'-r','LineWidth',2)
hold on
plot([freqcomp freqcomp],[0,1],'-k','LineWidth',2)
subplot(612);
title(sprintf('Z-2, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Coherence')
ylim([0 1])
plot([pper pper],[0 1],'-b','LineWidth',2)
hold on
plot([sper sper],[0 1],'-r','LineWidth',2)
hold on
plot([freqcomp freqcomp],[0,1],'-k','LineWidth',2)
subplot(613);
title(sprintf('P-Z, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Coherence')
ylim([0 1])
plot([pper pper],[0 1],'-b','LineWidth',2)
hold on
plot([sper sper],[0 1],'-r','LineWidth',2)
hold on
plot([freqcomp freqcomp],[0,1],'-k','LineWidth',2)
subplot(614);
title(sprintf('1-2, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Coherence')
ylim([0 1])
plot([pper pper],[0 1],'-b','LineWidth',2)
hold on
plot([sper sper],[0 1],'-r','LineWidth',2)
hold on
plot([freqcomp freqcomp],[0,1],'-k','LineWidth',2)
subplot(615);
title(sprintf('1-P, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Coherence')
ylim([0 1])
plot([pper pper],[0 1],'-b','LineWidth',2)
hold on
plot([sper sper],[0 1],'-r','LineWidth',2)
hold on
plot([freqcomp freqcomp],[0,1],'-k','LineWidth',2)
subplot(616);
title(sprintf('2-P, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Coherence')
ylim([0 1])
plot([pper pper],[0 1],'-b','LineWidth',2)
hold on
plot([sper sper],[0 1],'-r','LineWidth',2)
hold on
plot([freqcomp freqcomp],[0,1],'-k','LineWidth',2)
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);
filename=sprintf('%s/Coherence_%s',figoutpath,station);
print(gcf,'-dpng',filename)

figure(6)
subplot(611);
title(sprintf('Z-1, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Phase')
ylim([-90 90])
plot([pper pper],[-90 90],'-b','LineWidth',2)
hold on
plot([sper sper],[-90 90],'-r','LineWidth',2)
hold on
plot([freqcomp freqcomp],[-90 90],'-k','LineWidth',2)
hold on
subplot(612);
title(sprintf('Z-2, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Phase')
ylim([-90 90])
plot([pper pper],[-90 90],'-b','LineWidth',2)
hold on
plot([sper sper],[-90 90],'-r','LineWidth',2)
hold on
plot([freqcomp freqcomp],[-90 90],'-k','LineWidth',2)
hold on
subplot(613);
title(sprintf('P-Z, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Phase')
ylim([-90 90])
plot([pper pper],[-90 90],'-b','LineWidth',2)
hold on
plot([sper sper],[-90 90],'-r','LineWidth',2)
hold on
plot([freqcomp freqcomp],[-90 90],'-k','LineWidth',2)
hold on
subplot(614);
title(sprintf('1-2, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Phase')
ylim([-90 90])
plot([pper pper],[-90 90],'-b','LineWidth',2)
hold on
plot([sper sper],[-90 90],'-r','LineWidth',2)
hold on
plot([freqcomp freqcomp],[-90 90],'-k','LineWidth',2)
hold on
subplot(615);
title(sprintf('1-P, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Phase')
ylim([-90 90])
plot([pper pper],[-90 90],'-b','LineWidth',2)
hold on
plot([sper sper],[-90 90],'-r','LineWidth',2)
hold on
plot([freqcomp freqcomp],[-90 90],'-k','LineWidth',2)
hold on
subplot(616);
title(sprintf('2-P, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Phase')
ylim([-90 90])
plot([pper pper],[-90 90],'-b','LineWidth',2)
hold on
plot([sper sper],[-90 90],'-r','LineWidth',2)
hold on
plot([freqcomp freqcomp],[-90 90],'-k','LineWidth',2)
hold on
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);
filename=sprintf('%s/Phase_%s',figoutpath,station);
print(gcf,'-dpng',filename)

figure(5)
subplot(611);
title(sprintf('Z-1, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Admittance')
hold on
subplot(612);
title(sprintf('Z-2, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Admittance')
hold on
subplot(613);
title(sprintf('P-Z, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Admittance')
hold on
subplot(614);
title(sprintf('1-2, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Admittance')
hold on
subplot(615);
title(sprintf('1-P, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Admittance')
hold on
subplot(616);
title(sprintf('2-P, Station: %s',station))
xlabel('Frequency (Hz)')
ylabel('Admittance')
hold on
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);
filename=sprintf('%s/Admittance_%s',figoutpath,station);
print(gcf,'-dpng',filename)

% Calculate average spectra
czz_mean = mean(czz_all,1);
cpp_mean = mean(cpp_all,1);
c11_mean = mean(c11_all,1);
c22_mean = mean(c22_all,1);

czz_std = std(czz_all,1);
cpp_std = std(cpp_all,1);
c11_std = std(c11_all,1);
c22_std = std(c22_all,1);

czz_med = median(czz_all,1);
cpp_med = median(cpp_all,1);
c11_med = median(c11_all,1);
c22_med = median(c22_all,1);

% Calculate average coherence
coh1z_mean = mean(coh1z_all,1);
coh2z_mean = mean(coh2z_all,1);
cohpz_mean = mean(cohpz_all,1);
coh12_mean = mean(coh12_all,1);
coh1p_mean = mean(coh1p_all,1);
coh2p_mean = mean(coh2p_all,1);

coh1z_med = median(coh1z_all,1);
coh2z_med = median(coh2z_all,1);
cohpz_med = median(cohpz_all,1);
coh12_med = median(coh12_all,1);
coh1p_med = median(coh1p_all,1);
coh2p_med = median(coh2p_all,1);

coh1z_std = std(coh1z_all,1);
coh2z_std = std(coh2z_all,1);
cohpz_std = std(cohpz_all,1);
coh12_std = std(coh12_all,1);
coh1p_std = std(coh1p_all,1);
coh2p_std = std(coh2p_all,1);


filename = [outpath,station,'_spectraavg.mat'];
    save(filename,'station','f','czz_mean','cpp_mean','c11_mean','c22_mean','czz_std','cpp_std','c11_std','c22_std','czz_med','cpp_med','c11_med','c22_med',...
        'elev','freqcomp','coh1z_mean','coh2z_mean','cohpz_mean','coh12_mean','coh1p_mean','coh2p_mean','coh1z_med','coh2z_med','cohpz_med','coh12_med','coh1p_med','coh2p_med'...
        ,'coh1z_std','coh2z_std','cohpz_std','coh12_std','coh1p_std','coh2p_std');

