% script to go through data and correct for instrument response for  each 
%  of the station's channels. This script uses full response files from the
%  dataless seed file, and then the command line tool evalresp (from IRIS)
%  to compute the full response over a range of frequencies. 
clear all

addpath('matguts')

overwrite = true;

isfigure = 1; 

% path to instrument response data
respdir = '~/Work/CASCADIA/CAdb/RESP_NW_STA_CHA/';

% path to evalresp tool
evrespdir = '~/Work/Codes/evalresp-3.3.3/';

% antelope db details
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascattendb';

% path to top level of directory tree for data
datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash

% path to working directory
wdir = '/Users/zeilon/Documents/MATLAB/CASC_atten';
%% get to work

respfiles = dir(respdir); 

db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,elats,elons,edeps,evtimes,mags] = dbgetv(dbor,'orid','lat','lon','depth','time','ms');
norids = dbnrecs(dbor);
dbclose(db);

for ie = 196:196 % 1:norids % loop on orids
    fprintf('\n Orid %.0f %s \n\n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
    evdir = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
    ev_yyyy_jjj = epoch2str(evtimes(ie),'%Y %j');
    if ~exist([datadir,evdir,'_datinfo.mat'],'file'), fprintf('No data at all for this event\n'), continue, end
    load([datadir,evdir,'_datinfo.mat'])
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end

    for is = 1:length(datinfo) % loop on stas
        sta = datinfo(is).sta; % sta name
        fprintf('Station %-5s...',sta)
        if ~overwrite
            if datinfo(is).rmresp, fprintf(' done already\n'), continue, end % skip if already removed response
        end
        
        load([datadir,evdir,sta,'.mat']); % load sta data for this evt
        if ~overwrite
        	if data.rmresp, fprintf(' done already\n'), continue, end % skip if already removed response
        end
        
        chans = data.chans;
        chans_raw = data.raw.chans;
        
        % length of data, nyquist freq etc.
        T = data.tt(end)-data.tt(1);
        N = data.nsamps;
        
        % parse names of channels, if have been rotated and re-named
        ich = find(strcmp(chans.component,'H'));
        icn = find(strcmp(chans.component,'N'));
        ice = find(strcmp(chans.component,'E')); 
        icz = find(strcmp(chans.component,'Z'));
        if isempty(icn), icn = find(strcmp(chans.component,'1')); end
        if isempty(ice), ice = find(strcmp(chans.component,'2')); end
        ics = [ich icn ice icz];
        
        ichr = find(strcmp(chans_raw.component,'H'));
        icnr = find(strcmp(chans_raw.component,'N'));
        icer = find(strcmp(chans_raw.component,'E')); 
        iczr = find(strcmp(chans_raw.component,'Z'));
        if isempty(icnr), icnr = find(strcmp(chans_raw.component,'1')); end
        if isempty(icer), icer = find(strcmp(chans_raw.component,'2')); end
        icsr = [ichr icnr icer iczr];
        
        % loop on channels
        for ic = 1:length(chans.component)
            fprintf(' %s,',chans.name{ic})
            
            dat0 = data.dat(:,ic);
            dat0 = detrend(dat0);
            dat0 = flat_hanning_win(1:length(dat0),dat0,1,length(dat0),50);
            
            % find RESP file
            respfilestr = ['RESP.',data.network,'.',sta,'..',chans_raw.name{icsr == ics(ic)} ];
            respfile = dir([respdir,respfilestr]);
            if isempty(respfile), fprintf('NO RESP'), end
            % use evalresp to make AMP and PHASE files, whence complex response
            [ resp,faxis ] = mkresp_evalresp( respdir,data.network,sta,chans_raw.name{icsr == ics(ic)},T,N,ev_yyyy_jjj );
            w = 2*pi*faxis;
            
            if isfigure
                figure(33)
                clf
                set(gcf,'position',[360   514   900   400]);
                hold on
                subplot(1,2,1)
                set(gca,'fontsize',18)
                loglog(faxis,abs(resp),'rx');
                subplot(1,2,2)
                set(gca,'fontsize',18)
                semilogx(faxis,angle(resp),'rx');
            end
            
            % add high pass - to stop low freqs blowing up
            lo_corner = 0.005;  % in Hz
            npoles=5;
            
            lo_w=2*pi*lo_corner;
            hpfiltfrq=( ((w./lo_w).^(2*npoles))./(1+(w./lo_w).^(2*npoles)) );
            
            % combine to get transfer function
            norm_trans=hpfiltfrq(:)./resp(:);    % transfer function
            norm_trans(isnan(norm_trans)) = 0;

            % fourier transform, remove response, inverse transform
            fftdat0 = fft(dat0);
            fftdat0 = fftdat0(:).*norm_trans(:); % in f-space multiply (i.e. convolve time-domain function)
            dat0_cor = real(ifft(fftdat0));
            
            if isfigure
                figure(2); clf; hold on
                plot( data.dat(:,ic)./max(abs( data.dat(:,ic)) ) - 1,'b')
                plot( dat0          ./max(abs( dat0)           )    ,'r')
                plot( dat0_cor      ./max(abs( dat0_cor)       ) + 1,'g')
            end
            pause
            data.dat(:,ic) = dat0_cor;
        end % loop on chans
        
        %% log resp removal and save
        fprintf(' resp removed\n')
        data.rmresp = true;
        save([datadir,evdir,sta],'data')
        
        datinfo(is).rmresp = true;
        save([datadir,evdir,'_datinfo'],'datinfo')
            
    end % loop on stas

	%% sum up
    fprintf(' STA  CHAN  NEZ  resp  tilt  comp\n')
    for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f     %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).rmtilt,datinfo(is).rmcomp); end
            
end
