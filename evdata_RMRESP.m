% script to go through data and correct for instrumetn response for  each 
%  of the station's channels 
clear all
addpath('matguts')

overwrite = true;

isfigure = 0; 

% path to instrument response data
respdir = '~/Work/CASCADIA/CAdb/response/';

% antelope db details
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascattendb';

% path to top level of directory tree for data
datadir = '~/Work/CASCADIA/DATA/'; % needs final slash

%% get to work

respfiles = dir(respdir); 

db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,elats,elons,edeps,evtimes,mags] = dbgetv(dbor,'orid','lat','lon','depth','time','ms');
norids = dbnrecs(dbor);
dbclose(db);

for ie = 85:85 % 1:norids % loop on orids
    fprintf('\n Orid %.0f %s \n\n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
    evdir = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
    evday = epoch2str(evtimes(ie),'%Y%j');
%     stafiles = dir([datadir,evdir]); stafiles = stafiles(3:end);
    load([datadir,evdir,'_datinfo.mat'])
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end

    for is = 1:length(datinfo) % loop on stas
        sta = datinfo(is).sta; % sta name
        fprintf('Station %-5s...',sta)
        if datinfo(is).rmresp, fprintf(' done already\n'), continue, end % skip if already removed response
        
        db = dbopen([dbdir,dbnam],'r');
        dbsch = dblookup_table(db,'sitechan');
        dbschs = dbsubset(dbsch,sprintf('sta == "%s" && ondate <= %s && offdate >= %s',sta,evday,evday));
        [schans,respfilestrs] = dbgetv(dbschs,'chan','descrip');
        dbclose(db);
        if ~iscell(respfilestrs), respfilestrs = {respfilestrs}; end
        
        load([datadir,evdir,sta,'.mat']); % load sta data for this evt
        if data.rmresp, fprintf(' done already\n'), continue, end % skip if already removed response

        if strcmp(sta,'M02CO') || strcmp(sta,'M04CO'), sta = sta(1:4); end % rename problem stations
        
        chans = data.chans;
        chans_raw = data.raw.chans;
        
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
        
        for ic = 1:length(chans.component)
            fprintf(' %s,',chans.name{ic})
            % icsr == ics(ic) % Need this chicanery in case chnames changed
            respfilestri = respfilestrs{ strcmp( schans,chans_raw.name(icsr == ics(ic)) ) };
            respfile = dir([respdir,respfilestri,'*']);
            if isempty(respfile), fprintf('NO RESP'), end
            [zeros,poles,gain] = read_sac_pole_zero_haj([respdir,respfile.name]);
            
            %% RESPONSE REMOVAL, CRIBBED FROM JINGLE'S rm_resp SCRIPT
            % options
            lo_corner = 0.005;  % in Hz
            npoles=5;
            
            dat0 = data.dat(:,ic);
            dat0 = detrend(dat0);
            dat0 = flat_hanning_win(1:length(dat0),dat0,1,length(dat0),50);
            T = data.tt(end)-data.tt(1);
            N = data.nsamps;
            delta = 1./data.samprate;
            if mod(N,2)
                 faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/T);
            else
                 faxis = [0:N/2,-N/2+1:-1]*(1/T);
            end
            w = faxis.*2*pi;
            resp = ones(size(w));
            for ip = 1:length(poles)
                resp = resp./(1i*w - poles(ip));
            end
            for ip = 1:length(zeros)
                resp = resp.*(1i*w - zeros(ip));
            end
            resp = resp*gain;
            
            if isfigure
                figure(33)
                clf
                set(gcf,'position',[360   514   900   400]);
                hold on
                subplot(1,2,1)
                set(gca,'fontsize',18)
                semilogy(faxis,abs(resp),'rx');
                subplot(1,2,2)
                set(gca,'fontsize',18)
                plot(faxis,angle(resp),'rx');
            end
            
            lo_w=2*pi*lo_corner;
            hpfiltfrq=( ((w./lo_w).^(2*npoles))./(1+(w./lo_w).^(2*npoles)) );
            norm_trans=hpfiltfrq./resp;    % this is normalization transfer function
            norm_trans(isnan(norm_trans)) = 0;

            fftdat0 = fft(dat0);
            fftdat0 = fftdat0(:).*norm_trans(:);
            dat0_cor = real(ifft(fftdat0));
            
            if isfigure
                figure(2); clf; hold on
                plot( data.dat(:,ic)./max(abs(data.dat(:,ic))) - 1,'b')
                plot( dat0./max(abs(dat0)),'r')
                plot( dat0_cor./max(abs(dat0_cor)) + 1,'g')
            end
            
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
