% script to go through data and correct for instrumetn response for  each 
%  of the station's channels 
clear all
addpath('/Users/zeilon/Documents/MATLAB/helen_dataprocessing')
addpath('../matguts')

overwrite = false;

% day info to calculate tilt and compliance
noisewind = 43200;     % length of noise window before event time (sec).
T    = 6000;  % the length of each time window, in sec  
Nwin = 10;    % N of time windows into which to chop the noise 
if ~(T > noisewind/Nwin), error('Make T larger so there is overlap between windows'); end

fmin=0.005;fmax=0.02;
coh_min=0.8;
chorz_max=1;


isfigure = 0; 

% path to sensor tilt data
tiltdir = '~/Work/CASCADIA/CAdb/tilt_compliance/';

% antelope db details
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';

% path to top level of directory tree for data
datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash

%% get to work
[ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam );
[ nstas,stas,slats,slons,selevs,ondate,offdate,staname,statype ] = db_stadata( dbdir,dbnam );

for ie = 215:215 % 1:norids % loop on orids
    fprintf('\n Orid %.0f %s \n\n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
    evdir = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
    evday = epoch2str(evtimes(ie),'%Y%j');
    datinfofile = [datadir,evdir,'_datinfo'];
   
    % check files exist
    if ~exist([datinfofile,'.mat'],'file'), fprintf('No data at all for this event\n');continue, end
    load(datinfofile)
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end

    for is = 7:length(datinfo) % loop on stas
        sta = datinfo(is).sta; % sta name
        fprintf('Station %-5s...',sta)
        if datinfo(is).rmtilt, fprintf(' done already\n'), continue, end % skip if already removed tilt
        if strcmp(statype(strcmp(stas,sta)),'OBS') ~= 1, continue, end % skip if not an OBS
        
        load([datadir,evdir,sta,'.mat']); % load sta data for this evt
        if data.rmtilt, fprintf(' done already\n'), continue, end % skip if already removed tilt
        
        %% sort out channels        
        chans = data.chans;
        chans_raw = data.raw.chans;
        
        ich = find(strcmp(chans.component,'H'));
        icn = find(strcmp(chans.component,'N'));
        ice = find(strcmp(chans.component,'E')); 
        icz = find(strcmp(chans.component,'Z'));
        if isempty(icn), icn = find(strcmp(chans.component,'1')); end
        if isempty(ice), ice = find(strcmp(chans.component,'2')); end
        ics = [ich icn ice icz];
                
        if length(ics)~=4, continue, end % skip if don't have all 4 chans
        
        %% calc some basics
        samprate =data.samprate;
        dt = 1./data.samprate;
        npts = T*samprate;
        Ppower = nextpow2(npts);
        NFFT = 2^Ppower;
        npad0 = (NFFT-npts)/2;
        f = samprate/2*linspace(0,1,NFFT/2+1);
        win_t0 = round(linspace(1,noisewind-T+1,Nwin));
        win_pt0 = (win_t0-1)*samprate+1;

        %% scales for seismometer and pressure chans...?
        scalez = 1;
        scalep = 1;
        
        %% Prep for spectra
        % ----- open zeros array for spectrum and conj(spectrum) for 4 channels -----
        spectrum_Z=zeros(Nwin,length(f));
        spectrum_H1=zeros(Nwin,length(f));
        spectrum_H2=zeros(Nwin,length(f));
        spectrum_P=zeros(Nwin,length(f));
        cspectrum_Z=zeros(Nwin,length(f));
        cspectrum_H1=zeros(Nwin,length(f));
        cspectrum_H2=zeros(Nwin,length(f));
        cspectrum_P=zeros(Nwin,length(f));
        
        
        for iwin = 1:Nwin % Window index
            j1 = win_pt0(iwin);
            j2 = j1+npts-1;

            datZ  = data.dat(j1:j2,icz);
            datZ(isnan(datZ)) = 0;
            datZ  = datZ.*flat_hanning_taper(data.tt(j1:j2),0.4*dt*npts);
            datZ  = detrend(datZ,0);
            %amp_Z  = detrend(amp_Z,'linear');
            %amp_Z  = cos_taper(amp_Z);
            %amp_Z  = padarray(amp_Z,[0 npad0],'both');
            spectrum = fft(datZ,NFFT); 
            spectrum_Z(iwin,1:length(f)) = spectrum(1:NFFT/2+1)/scalez;
            cspectrum_Z(iwin,1:length(f)) = conj(spectrum(1:NFFT/2+1))/scalez;
%             [spect_Z,FZ,TZ] = spectrogram(datZ,flat_hanning(data.tt(j1:j2),0.4*dt*npts),0,npts,samprate,'yaxis');

            %loglog(f,2*abs(spectrum(1:NFFT/2+1))) 
            %hold on;
            %title('Single-Sided Amplitude Spectrum of y(t)')

            datN = data.dat(j1:j2,icn);
            datN(isnan(datN)) = 0;
            datN = datN.*flat_hanning_taper(data.tt(j1:j2),0.4*dt*npts);
            datN = detrend(datN,0);
            %amp_H1 = detrend(amp_H1,'linear');
            %amp_H1 = cos_taper(amp_H1);
            %amp_H1  = padarray(amp_H1,[0 npad0],'both');
            spectrum = fft(datN,NFFT);
            spectrum_H1(iwin,1:length(f)) = spectrum(1:NFFT/2+1)/scalez;
            cspectrum_H1(iwin,1:length(f)) = conj(spectrum(1:NFFT/2+1))/scalez;
			
%             %% --- quick plot check abs(spectrum)^2 == spectrum*cspectrum 
% 			loglog(f,(abs(spectrum(1:NFFT/2+1))).^2) 
% 			pause
% 			hold on;
% 			loglog(f,spectrum(1:NFFT/2+1) .* conj(spectrum(1:NFFT/2+1)),'r-')
% 			% --- end quick plot
            
            datE = data.dat(j1:j2,ice);
            datE(isnan(datE)) = 0;
            datE = datE.*flat_hanning_taper(data.tt(j1:j2),0.4*dt*npts);
            datE = detrend(datE,0);
            %amp_H2 = detrend(amp_H2,'linear');
            %amp_H2 = cos_taper(amp_H2);
            %amp_H2  = padarray(amp_H2,[0 npad0],'both');
            spectrum = fft(datE,NFFT);
            spectrum_H2(iwin,1:length(f)) = spectrum(1:NFFT/2+1)/scalez;
            cspectrum_H2(iwin,1:length(f)) = conj(spectrum(1:NFFT/2+1))/scalez;
            %loglog(f,2*abs(spectrum(1:NFFT/2+1))) 
            %hold on;

            datP  = data.dat(j1:j2,ich);
            datP(isnan(datP)) = 0;
            datP  = datP.*flat_hanning_taper(data.tt(j1:j2),0.4*dt*npts);
            datP  = detrend(datP,0);
            %amp_P  = detrend(amp_P,'linear');
            %amp_P  = cos_taper(amp_P);
            %amp_P  = padarray(amp_P,[0 npad0],'both');
            spectrum = fft(datP,NFFT);
            spectrum_P(iwin,1:length(f)) = spectrum(1:NFFT/2+1)/scalep; 
            cspectrum_P(iwin,1:length(f)) = conj(spectrum(1:NFFT/2+1))/scalep;
            %loglog(f,2*abs(spectrum(1:NFFT/2+1))) 
            %hold on;
                
                
        end % for iwin
            
        %% compute auto and cross spectra
        % ---------- get auto-spectra g11,g22,gPP, gZZ 
        %            and cross-spectra g1Z,	g2Z,gPZ, g12, g1P, g2P 
        %            similar to eq(6) in Crawford and Webb, BSSA2000 ----------
        g11 = spectrum_H1.*cspectrum_H1;
        g11 = 2/Nwin/(NFFT*dt)*sum(g11,1); % average across all winds

        g22 = spectrum_H2.*cspectrum_H2;
        g22 = 2/Nwin/(NFFT*dt)*sum(g22,1); % average across all winds

        gPP = spectrum_P.*cspectrum_P;
        gPP = 2/Nwin/(NFFT*dt)*sum(gPP,1); % average across all winds

        gZZ = spectrum_Z.*cspectrum_Z;
        gZZ = 2/Nwin/(NFFT*dt)*sum(gZZ,1); % average across all winds

        g1Z = cspectrum_H1 .* spectrum_Z;
        g1Z = 2/Nwin/(NFFT*dt)*sum(g1Z,1); % average across all winds

        g2Z = cspectrum_H2 .* spectrum_Z;
        g2Z = 2/Nwin/(NFFT*dt)*sum(g2Z,1); % average across all winds  

        gPZ = cspectrum_P .* spectrum_Z;
        gPZ = 2/Nwin/(NFFT*dt)*sum(gPZ,1); % average across all winds    

        g12 = cspectrum_H1 .* spectrum_H2;
        g12 = 2/Nwin/(NFFT*dt)*sum(g12,1); % average across all winds

        g1P = cspectrum_H1 .* spectrum_P;
        g1P = 2/Nwin/(NFFT*dt)*sum(g1P,1); % average across all winds

        g2P = cspectrum_H2 .* spectrum_P;
        g2P = 2/Nwin/(NFFT*dt)*sum(g2P,1); % average across all winds
        
        %% plot coherences before removal
            
        cohhj_z1=abs(g1Z).^2./(g11.*gZZ);
        cohhj_z2=abs(g2Z).^2./(g22.*gZZ);
        cohhj_zp=abs(gPZ).^2./(gPP.*gZZ);
        cohhj_12=abs(g12).^2./(g11.*g22);
        cohhj_1p=abs(g1P).^2./(g11.*gPP);
        cohhj_2p=abs(g2P).^2./(g22.*gPP);
        figure(88)
        semilogx(f,abs(cohhj_z1),'r'); hold on
        semilogx(f,abs(cohhj_z2),'b'); hold on
        semilogx(f,abs(cohhj_zp),'k')
        legend('Z1','Z2','ZP');
            
%             compliance_PS = [figoutpath,char(eventids(ie)),'_',char(ioridEVT(ie)),'_',staname,'_coherence1.ps'];
%         print('-dpsc2',compliance_PS);
            
            figure(89)
            semilogx(f,abs(cohhj_12),'r'); hold on
            semilogx(f,abs(cohhj_1p),'b'); hold on
            semilogx(f,abs(cohhj_2p),'k')
            legend('12','1p','2p');
%             
%             compliance_PS = [figoutpath,char(eventids(ie)),'_',char(ioridEVT(ie)),'_',staname,'_coherence2.ps'];
%         print('-dpsc2',compliance_PS);
        
        
        % ----------get the ZZ auto-spectrum after removing channel 1 gZZ_1, 
        %           then 2 gZZ_12, and then P, gZZ_12P
        %     plus  coherence between P and Z after removing effects of 1 and 2, gamPZ_12
        %     plus  ransfer function, l1Z,l12,l1P,l2P_1,l2Z_1,lPZ_12 
        %                             l1Z == conj(AZ1) in eq(8) in Crawford and Webb, BSSA2000 ----------
%         [gZZ_12P,gZZ_12,gZZ_1,gamPZ_12,l1Z,l12,l1P,l2P_1,l2Z_1,lPZ_12]=multicoher_p(g11,g22,gPP,gZZ,g1Z,g2Z,gPZ,g12,g1P,g2P,f);
        [gZZ_12P,gZZ_12,gZZ_1,gamPZ_12,l1Z,l12,l1P,l2P_1,l2Z_1,lPZ_12]=multicoher(g11,g22,gPP,gZZ,g1Z,g2Z,gPZ,g12,g1P,g2P,f);
        w=2*pi*f;
        gZZ_12P=gZZ_12P./(1i*w);

        % ---------- plot PS files for auto-Spectrum for gZZ, gZZ_12, and gZZ_12P ---------
        figure(102)
        clf
        loglog(f,abs(gZZ),'k-');hold on;
        loglog(f,abs(gZZ_12),'r-');hold on;
        loglog(f,abs(gZZ_12P),'b-');hold on;
        legend('gzz','gZZ-12','gZZ-12P')
        
%         title([char(eventids(ie)),'_',char(ioridEVT(ie)),'-',staname,'-auto-Spectrum for ZZ Nwin ',num2str(Nwin),' in a day, in ', num2str(nday),' days']);
        xlim([0.001 20]);
        return
        
%         AutoSpectrumZZ_PS = [figoutpath,char(eventids(ie)),'_',char(ioridEVT(ie)),'_',staname,'_AutoSpectrumZZ.ps'];
%         print('-dpsc2',AutoSpectrumZZ_PS);

        
        figure(103)
        clf
        loglog(f,abs(g11),'b-');hold on;
        loglog(f,abs(g22),'r-');hold on;
        loglog(f,abs(gZZ),'k-');hold on;
        legend('g11','g22','gzz')
        
%         title([char(eventids(ie)),'_',char(ioridEVT(ie)),'-',staname,'-auto-Spectrum for ZZ Nwin ',num2str(Nwin),' in a day, in ', num2str(nday),' days']);
        xlim([0.001 20]);
        
%         AutoSpectrumZZ_PS = [figoutpath,char(eventids(ie)),'_',char(ioridEVT(ie)),'_',staname,'_AutoSpectrum1122.ps'];
%         print('-dpsc2',AutoSpectrumZZ_PS);
        
        figure(104)
        clf
        loglog(f,abs(g11),'b-');hold on;
        loglog(f,abs(g22),'r-');hold on;
        loglog(f,abs(gZZ),'k-');hold on;
        legend('g11','g22','gzz')
        
%         title([char(eventids(ie)),'_',char(ioridEVT(ie)),'-',staname,'-auto-Spectrum for ZZ Nwin ',num2str(Nwin),' in a day, in ', num2str(nday),' days']);
        xlim([0.05 1]);
        
        AutoSpectrumZZ_PS = [figoutpath,char(eventids(ie)),'_',char(ioridEVT(ie)),'_',staname,'_AutoSpectrum1122_zoom.ps'];
        print('-dpsc2',AutoSpectrumZZ_PS);
        
        figure(105) 
        subplot(2,1,1)
		loglog(f,smooth(abs(lPZ_12),100),'b-');
        xlim([0.005 25]);
%         title([char(eventids(ie)),'_',char(ioridEVT(ie)),'-',staname,'-compliance: transfer function between P and Z after removing the effects of channels 1 and 2 '] )
        subplot(2,1,2)
	semilogx(f,smooth(gamPZ_12,100));	
        title('Coherence between P and Z');
	xlim([0.005 25]);
        
        
        return
        %% log resp removal and save
        fprintf(' resp removed\n')
        data.rmresp = true;
        save([datadir,evdir,sta],'data')
        
        datinfo(is).rmresp = true;
        save(datinfofile,'datinfo')
                    
    end % loop on stas

	%% sum up
    fprintf(' STA  CHAN  NEZ  resp  tilt  comp\n')
    for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f     %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).rmtilt,datinfo(is).rmcomp); end
            
end
