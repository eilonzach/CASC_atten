% Script to compute spectral characteristics of noise before each event.
% This script goes through data event by event, then station by station.
% For OBS stations, it pulls out the 12 hours of noise prior to each event
% and computes spectra on a set of Nwin windows of noise, saving the
% outputs into structures within the data directories
clear all

overwrite = true;

resamprate = 5; % NEW SAMPLE RATE TO DOWNSAMP TO

% different methods to compute transfer function - from Helen
% 1 = TILT AND COMPLIANCE, CRAWFORD AND WEBB, 2000, spahr's code
% 2 = COMPLIANCE ONLY, same nomenclature but only removing pressure component
% 3 = COMPLIANCE AND THEN TILT, individual component method
% (4) = TILT AND COMPLIANCE - BELL ET AL., 2014 reorient in maximum coherence direction
method_trans_func = 1; 

fmin=0.005;fmax=0.02;
coh_min=0.8;
chorz_max=1;

pper = 1/16;
sper = 1/8;



ifplot = 1; 
ifsavefigs = 1;

%% paths
cd('/Users/zeilon/Documents/MATLAB/CASC_atten/')
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

for ie = 215:215 % 1:norids % loop on orids
    fprintf('\n Orid %.0f %s \n\n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
    evdir = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
    evday = epoch2str(evtimes(ie),'%Y%j');
    datinfofile = [datadir,evdir,'_datinfo'];
           
    if any((evtimes(ie)-evtimes)>0 & (evtimes(ie)-evtimes) < 24*60*60)
        fprintf('Another event within prev 24 hrs... skipping\n')
        continue
    end

    % check files exist
    if exist([datinfofile,'.mat'],'file')~=2, fprintf('No data at all for this event\n');continue, end
    load(datinfofile)
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end

    for is = 1:length(datinfo) % loop on stas
        sta = datinfo(is).sta; % sta name
        
        %% load data, check for previous activity
        if strcmp(statype(strcmp(stas,sta)),'OBS') ~= 1, continue, end % skip if not an OBS
        
        specfile = [datadir,evdir,sta,'_spectra.mat'];
        datfile = [datadir,evdir,sta,'.mat'];
        
        % check files exist, otherwise skip
        if exist(specfile,'file')~=2, fprintf('no spectra file\n'), continue, end
        if exist(datfile,'file')~=2, fprintf('no data file\n'),  continue, end
        % conditions before loading - don't unecessarily load
        if datinfo(is).rmtilt && datinfo(is).rmtilt && ~overwrite
            fprintf(' done already\n'), continue
        end

        % OKAY - load data
        fprintf('Station %.0f, %-5s...',is,sta)
        load(datfile)
        load(specfile)
        
        % second pass on conditions
        if data.rmtilt && data.rmcomp && ~overwrite, fprintf(' done already\n'), continue, end % skip if already removed tilt + comp
        
        %% calc some basics
        T = (win_pt_end(1)-win_pt_start(1)+1)*dt;
        Nwin = length(win_pt_start);
        samprate = 1./dt;
        freqcomp = sqrt(9.8/(-2*pi*selevs(is)));

        
        %% COMPUTE TRANSFER FUNCTION
        
        switch method_trans_func
            case 1
                [czz_12p,czz_12,czz_1,cohpz_12,l1z,l12,l1p,l2p_1,l2z_1,lpz_12,coh2z_1,cohpz_1,coh2p_1]...
                    = multicoher(c11_stack,c22_stack,cpp_stack,czz_stack,c1z_stack,c2z_stack,cpz_stack,c12_stack,c1p_stack,c2p_stack,f);
                transfer_fun = struct('sta',sta,'f',f,'freqcomp',freqcomp,'NFFT',NFFT,...
                                    'czz_12p',czz_12p,'czz_12',czz_12,'czz_1',czz_1,...
                                    'l1z',l1z,'l12',l12,'l1p',l1p,'l2p_1',l2p_1,'l2z_1',l2z_1,'lpz_12',lpz_12,...
                                    'cohpz_12',cohpz_12,'coh2z_1',coh2z_1,'cohpz_1',cohpz_1,'coh2p_1',coh2p_1);
            case 2
                [czz_p,lpz,lp1,lp2,coh1z_p,coh2z_p,coh12_p,c11_p,c22_p,c1z_p,c2z_p,c12_p]...
                    = compcohere(c11_stack,c22_stack,cpp_stack,czz_stack,c1z_stack,c2z_stack,cpz_stack,c12_stack,c1p_stack,c2p_stack,f);
                transfer_fun = struct('sta',sta,'f',f,'freqcomp',freqcomp,'NFFT',NFFT,...
                                    'czz_p',czz_p,'lpz',lpz,'lp1',lp1,'lp2',lp2,...
                                    'coh1z_p',coh1z_p,'coh2z_p',coh2z_p,'coh12_p',coh12_p,...
                                    'c11_p',c11_p,'c22_p',c22_p,'c1z_p',c1z_p,'c2z_p',c2z_p,'c12_p',c12_p);
            case 3
                [czz_p21,czz_p2,czz_p,coh1z_p2,lpz,lp2,lp1,l21_p,l2z_p,l1z_p2,coh2z_p,coh1z_p,coh21_p]...
                    = multicoher(cpp_stack,c22_stack,c11_stack,czz_stack,cpz_stack,c2z_stack,c1z_stack,c2p_stack,c1p_stack,c12_stack,f);
                transfer_fun = struct('sta',sta,'f',f,'freqcomp',freqcomp,'NFFT',NFFT,...
                                    'czz_p21',czz_p21,'czz_p2',czz_p2,'czz_p',czz_p,...
                                    'lpz',lpz,'lp1',lp1,'lp2',lp2,'l21_p',l21_p,'l2z_p',l2z_p,'l1z_p2',l1z_p2,...
                                    'coh1z_p2',coh1z_p2,'coh2z_p',coh2z_p,'coh1z_p',coh1z_p,'coh21_p',coh21_p);
            case 4
                fprintf('Sorry, this method not yet ready for primetime!\n'), return
        end
        if ifplot
            plot_trans_func(method_trans_func,transfer_fun) %< option in here to save, can switch on
        end

        
        %% CORRECT THE SEISMOGRAMS
        
        if pper < freqcomp
            lp=pper; %low pass filter transfer function at primary microseism
        else 
            lp=freqcomp;
        end        
        
        switch method_trans_func
            case 1
                % Filter transfer function to avoid microseism band
                l1z=gauss2_freq(l1z',f',lp,NFFT);
                l12=gauss2_freq(l12',f',lp,NFFT);
                l1p=gauss2_freq(l1p',f',lp,NFFT);
                l2p_1=gauss2_freq(l2p_1',f',lp,NFFT);
                l2z_1=gauss2_freq(l2z_1',f',lp,NFFT);
                lpz_12=gauss2_freq(lpz_12',f',lp,NFFT);
            case 2
                % Filter transfer function to avoid microseism band
                [lpz]=gauss2_freq(lpz',f',lp,NFFT);
            case 3
                % rename things
                l1z = lpz;
                l12 = lp2;
                l1p = lp1;
                l2p_1 = l21_p;
                l2z_1 = l2z_p;
                lpz_12 = l1z_p2;
                
                % Filter transfer function to avoid microseism band
                l1z=gauss2_freq(l1z',f',lp,NFFT);
                l12=gauss2_freq(l12',f',lp,NFFT);
                l1p=gauss2_freq(l1p',f',lp,NFFT);
                l2p_1=gauss2_freq(l2p_1',f',lp,NFFT);
                l2z_1=gauss2_freq(l2z_1',f',lp,NFFT);
                lpz_12=gauss2_freq(lpz_12',f',lp,NFFT);
        end
        
        %% %%%%%%%%% PROCESS EARTHQUAKE DATA
        
        % sort out channels        
        chans = data.chans;        
        ich = find(strcmp(chans.component,'H'));
        icn = find(strcmp(chans.component,'N'));
        ice = find(strcmp(chans.component,'E')); 
        icz = find(strcmp(chans.component,'Z'));
        if isempty(icn), icn = find(strcmp(chans.component,'1')); end
        if isempty(ice), ice = find(strcmp(chans.component,'2')); end
        ics = [ich icn ice icz];
        if length(ics)~=4  % skip if don't have all 4 chans
            fprintf('only %.0f chans... skipping\n',length(ics));  
            continue, 
        end 
        
        % grab data
        Zdat = data.dat(:,icz);
        H1dat = data.dat(:,icn);
        H2dat = data.dat(:,ice);
        Pdat = data.dat(:,ich);
        
        % NEED TO USE THE RAW PRESSURE IF THERE IS A LAMONT APG
        if strcmp(which_OBS(sta),'LDEO')
        	ichr = find(strcmp(raw.chans.component,'H'));
            Pdat = data.raw.dat(:,ichr);
            lo_corner = 0.005;  % in Hz
            npoles=5;
            lo_w=2*pi*lo_corner;

            N = length(Pdat);
            delta = dtP;
            Tr = N*delta;

            if mod(N,2)
                faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/Tr);
            else
                faxis = [0:N/2,-N/2+1:-1]*(1/Tr);
            end
            w = faxis.*2*pi;

            hpfiltfrq=( ((w./lo_w).^(2*npoles))./(1+(w./lo_w).^(2*npoles)) );
            norm_trans=hpfiltfrq;    % this is normalization transfer function
            norm_trans(isnan(norm_trans)) = 0;

            fftdata = fft(Pdat);
            fftdata = fftdata(:).*norm_trans(:);
            Pdat = real(ifft(fftdata));
        end
        
        %%Plot event data to check all components are working during event time
        fn = 1/2/dt;
        T1 = 20; T2= 100;
        [b,a]=butter(2,[1/fn/T2,1/fn/T1]);

        Z_filt  = filtfilt(b,a,Zdat);
        H1_filt  = filtfilt(b,a,H1dat);
        H2_filt  = filtfilt(b,a,H2dat);
        P_filt  = filtfilt(b,a,Pdat);

        figure(101)
        clf
        subplot(411)
        plot(taxisZ,Z_filt,'-k');
        hold on
        title(sprintf('%s Z',eventid));
        subplot(412)
        plot(taxisZ,H1_filt,'-k');
        hold on
        title(sprintf('%s H1',eventid));
        subplot(413)
        plot(taxisZ,H2_filt,'-k');
        hold on
        title(sprintf('%s H2',eventid));
        subplot(414)
        plot(taxisZ,P_filt,'-k');
        hold on
        title(sprintf('%s P',eventid));

        set(gcf,'PaperPositionMode','manual');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperOrientation','portrait');
        set(gcf,'PaperPosition',[.05 .05 8 10.5]);
    
        
                
        return
                
        %% SAVE
%         save(ofile,'sta','evdir','f','tt','dt','NFFT','win_pt_start','win_pt_end',...
%             'c11_stack','c22_stack','czz_stack','cpp_stack',...
%             'c12_stack','c1z_stack','c2z_stack','c1p_stack','c2p_stack','cpz_stack',...
%             'C12_stack','C1z_stack','C2z_stack','C1p_stack','C2p_stack','Cpz_stack',...
%             'Q12_stack','Q1z_stack','Q2z_stack','Q1p_stack','Q2p_stack','Qpz_stack',...
%             'coh12_stack','coh1z_stack','coh2z_stack','coh1p_stack','coh2p_stack','cohpz_stack',...
%             'ph12_stack','ph1z_stack','ph2z_stack','ph1p_stack','ph2p_stack','phpz_stack',...
%             'ad12_stack','ad1z_stack','ad2z_stack','ad1p_stack','ad2p_stack','adpz_stack',...
%             'spectrum_P','spectrum_Z','spectrum_H1','spectrum_H2',...
%             'cspectrum_P','cspectrum_Z','cspectrum_H1','cspectrum_H2');
%         fprintf('spectral details saved.\n')
%           
%         
%         %% plot noiseless data
%         if ifplot
%         plot_1evt1sta_spectra(ofile,ifsavefigs,[dbdir,'spectra/',sta])
% %         pause
%         end
%         %% log spectra calc and save
% 
%         datinfo(is).spectra = true;
%         save(datinfofile,'datinfo')
%                     
    end % loop on stas

	%% sum up
    fprintf(' STA  CHAN  NEZ  resp  spectra  tilt  comp\n')
    for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f      %1.0f       %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).spectra,datinfo(is).rmtilt,datinfo(is).rmcomp); end
            
end
