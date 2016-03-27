% Script to compute spectral characteristics of noise before each event.
% This script goes through data event by event, then station by station.
% For OBS stations, it pulls out the 12 hours of noise prior to each event
% and computes spectra on a set of Nwin windows of noise, saving the
% outputs into structures within the data directories
clear all

%% parms
overwrite = true;

ifplot = 1; 
ifsavefigs = 1;

%% paths
cd('/Users/zeilon/Documents/MATLAB/CASC_atten')
addpath('/Users/zeilon/Documents/MATLAB/helen_dataprocessing')
addpath('matguts')
% path to sensor spectra data
specdir = '~/Work/CASCADIA/CAdb/spectra/';
% antelope db details
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% path to top level of directory tree for data
datadir = '/Volumes/DATA_mini/CASCADIA/DATA/'; % needs final slash



%% get to work
[ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam ); % load events data
[ nstas,stas,slats,slons,selevs,ondate,offdate,staname,statype ] = db_stadata( dbdir,dbnam ); %load stations data

for is = 121:nstas % loop on stas
    sta = stas{is}; % sta name
   	if strcmp(statype(is),'OBS') ~= 1, continue, end % skip if not an OBS
	
    fprintf('Station %.0f, %-5s... ',is,sta)
    ofile = [specdir,sta,'/',sta,'_spectraavg.mat'];
    
	if exist(ofile,'file')==2
        fprintf('avg spectra file exists ')
        if overwrite, fprintf('- overwriting'), end
        if ~overwrite, fprintf('- skipping\n'), continue, end
    end
    
    %% start loop on orids
    kk = 0;
    for ie = 1:norids % 1:norids % loop on orids
        evdir = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];

        specfile = [datadir,evdir,sta,'_spectra.mat'];
        % check files exist
        if exist(specfile,'file')~=2, continue, end

        fprintf('\n Orid %.0f %s',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
        
        %% load specfile
        load(specfile);
        kk = kk+1;
        freqcomp = sqrt(9.8/(-2*pi*selevs(is)));

%         if check_spec ==1;
%         in1 = input('Would you like to use data from this day (n/y)? ','s');
%         if in1 == 'n';
%             in2 = input('Are you sure you want to delete this day (y)? ','s');
%             if in2 == 'y'
%                 delete(fullfile(inpath,spectra_filenames(ie).name));
%                 continue
%             end
%         end
%         end
        

        maxpowz(kk) = max(czz_stack)*10;
        maxpow1(kk) = max(c11_stack)*10;
        maxpow2(kk) = max(c22_stack)*10;
        maxpowp(kk) = max(cpp_stack)*10;

        minpowz(kk) = min(czz_stack)/10;
        minpow1(kk) = min(c11_stack)/10;
        minpow2(kk) = min(c22_stack)/10;
        minpowp(kk) = min(cpp_stack)/10;

        czz_all(kk,:) = czz_stack;
        cpp_all(kk,:) = cpp_stack;
        c11_all(kk,:) = c11_stack;
        c22_all(kk,:) = c22_stack;

        coh1z_all(kk,:) = coh1z_stack;
        coh2z_all(kk,:) = coh2z_stack;
        cohpz_all(kk,:) = cohpz_stack;
        coh12_all(kk,:) = coh12_stack;
        coh1p_all(kk,:) = coh1p_stack;
        coh2p_all(kk,:) = coh2p_stack;

    end % loop on events
    
    if kk == 0 % this station has no events
        fprintf('no evt spectra... skipping\n')
        continue
    end
    
    fprintf('\n computing averages... ')
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

    
    save(ofile,'sta','f','freqcomp',...
        'czz_mean','cpp_mean','c11_mean','c22_mean',...
        'czz_std', 'cpp_std', 'c11_std', 'c22_std',...
        'czz_med', 'cpp_med', 'c11_med', 'c22_med',...
        'coh1z_mean','coh2z_mean','cohpz_mean','coh12_mean','coh1p_mean','coh2p_mean',...
        'coh1z_med', 'coh2z_med', 'cohpz_med', 'coh12_med', 'coh1p_med', 'coh2p_med',...
        'coh1z_std', 'coh2z_std', 'cohpz_std', 'coh12_std', 'coh1p_std', 'coh2p_std');
    fprintf('mean spectra saved\n')

    %% plot
    if ifplot
        plot_stav_spectra
    end
    
end % loop on stas
