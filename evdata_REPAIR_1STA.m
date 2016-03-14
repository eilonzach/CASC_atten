% script to go through data and reset to raw traces and channels=
clear all
cd ~/Documents/MATLAB/CASC_atten/
addpath('matguts')

station = 'J44C';
phases = {'P','S','SKS','PKS'};
orientation = 300;


ifsave = true;

% antelope db details
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% path to top level of directory tree for data
datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash

%% get to work
db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,evtimes] = dbgetv(dbor,'orid','time');
norids = dbnrecs(dbor);

dbsch = dblookup_table(db,'sitechan');
dbsen = dblookup_table(db,'sensor');
dbins = dblookup_table(db,'instrument');        
%         dbschs = dbsubset(dbsch,sprintf('sta == "%s" && ondate <= %s && offdate >= %s',sta,evday,evday));
dbschs = dbsubset(dbsch,sprintf('sta == "%s"',station));
dbj1 = dbjoin(dbschs,dbsen);
dbj2 = dbjoin(dbj1,dbins);
[schans,respdirs,respfiles] = dbgetv(dbj2,'chan','dir','dfile');
dbclose(db);

for ie = 224:norids % 1:norids % loop on orids
    fprintf('Orid %.0f %s \n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
    evdir    = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
    datinfofile = [datadir,evdir,'_datinfo'];
    load(datinfofile)
    if isempty(intersect({datinfo.sta}',station)), continue; end % sta has no data for this evt
    is = find(strcmp({datinfo.sta}',station));
    
    load([datadir,evdir,station,'.mat']); stadat = data;% load sta data for this evt
    chans = stadat.chans;
    raw = stadat.raw;
    
    %% do any corrections
    
    
    %% redo rotation
    ich = find(strcmp(chans.component,'H'));
    icn = find(strcmp(chans.component,'N'));if isempty(icn), icn = find(strcmp(chans.component,'1')); end
    ice = find(strcmp(chans.component,'E'));if isempty(ice), ice = find(strcmp(chans.component,'2')); end 
    icz = find(strcmp(chans.component,'Z'));

    ics = [ich icn ice icz];

    ichr = find(strcmp(raw.chans.component,'H'));
    icnr = find(strcmp(raw.chans.component,'N'));if isempty(icnr), icnr = find(strcmp(raw.chans.component,'1')); end
    icer = find(strcmp(raw.chans.component,'E'));if isempty(icer), icer = find(strcmp(raw.chans.component,'2')); end
    iczr = find(strcmp(raw.chans.component,'Z'));
    
    ic_new = [ichr icnr icer iczr];

    newdat = zeros(size(stadat.dat));
    sor = sind(orientation);
    cor = cosd(orientation);

    newdat(:,icn) = cor*raw.dat(:,icn) - sor*raw.dat(:,ice);    % NORTH
    newdat(:,ice) = sor*raw.dat(:,icn) + cor*raw.dat(:,ice);    % EAST
    newdat(:,ich) = raw.dat(:,ich);
    newdat(:,icz) = raw.dat(:,icz);

    newdat = newdat(:,ic_new);
    
    %% redo response removal
    for ic = 1:length(chans.component)
        fprintf(' %s,',chans.name{ic})
        % icsr == ics(ic) % Need this chicanery in case chnames changed
        iresp = strcmp( schans,raw.chans.name(ic_new == ics(ic)) ) ;
%         iresp = strcmp( schans,chans.name(ic) ) ;
        respfile = dir([respdirs{iresp},respfiles{iresp},'*']);
        if isempty(respfile), fprintf('NO RESP'), end
        [zz,pp,gain] = read_sac_pole_zero([respdirs{iresp},respfiles{iresp}]);
        %% RESPONSE REMOVAL, CRIBBED FROM JINGLE'S rm_resp SCRIPT
        lo_corner = 0.005;  % in Hz
        npoles=5;

        dat0 = newdat(:,ic);
        dat0 = detrend(dat0);
        dat0 = flat_hanning_win(1:length(dat0),dat0,1,length(dat0),50);

        T = stadat.tt(end)-stadat.tt(1);
        N = stadat.nsamps;
        delta = 1./stadat.samprate;
        if mod(N,2)
             faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/T);
        else
             faxis = [0:N/2,-N/2+1:-1]*(1/T);
        end
        w = faxis.*2*pi;
        resp = ones(size(w));
        for ip = 1:length(pp)
            resp = resp./(1i*w - pp(ip));
        end
        for ip = 1:length(zz)
            resp = resp.*(1i*w - zz(ip));
        end
        resp = resp*abs(gain); % made abs to get rid of neg gain issues!

        lo_w=2*pi*lo_corner;
        hpfiltfrq=( ((w./lo_w).^(2*npoles))./(1+(w./lo_w).^(2*npoles)) );
        norm_trans=hpfiltfrq./resp;    % transfer function
        norm_trans(isnan(norm_trans)) = 0;

        fftdat0 = fft(dat0);
        fftdat0 = fftdat0(:).*norm_trans(:); % in f-space multiply (i.e. convolve time-domain function)
        dat0_cor = real(ifft(fftdat0));

        newdat(:,ic) = dat0_cor;
    end % loop on chans
        
        
    %% fix eqar structure
    stadat.dat    = newdat;
    stadat.NEZ    = true;
    stadat.rmresp = true;
    stadat.rmtilt = false;
    stadat.rmcomp = false;
    
    [stadat.phases.xtime]   = deal([]);
    [stadat.phases.xartime] = deal([]);
    [stadat.phases.xacor]   = deal([]);
    
    if ifsave
        data = stadat;
        save([datadir,evdir,station],'data')
    end
    
    %% fix datinfo
    datinfo(is).NEZ    = true;
    datinfo(is).rmresp = true;
    datinfo(is).rmtilt = false;
    datinfo(is).rmcomp = false;
    if ifsave
        save(datinfofile,'datinfo')
    end
    %% fix phase-specific structures
    fprintf('Phases: ')
    for ip = 1:length(phases)
        phase = phases{ip};
        fprintf('%s... ',phase)
        
        phdatinfofile = [datadir,evdir,'_datinfo_',phase,'.mat'];
        if isempty(dir(phdatinfofile)), continue, end
        load(phdatinfofile)
        is = find(strcmp({datinfo.sta}',station));
             
        eqarchans = {'Z','R','T'};
        for ic = 1:length(eqarchans)
            fprintf('%s',eqarchans{ic})
            eqarfile = [datadir,evdir,'_EQAR_',phase,'_',eqarchans{ic},'.mat'];
            load(eqarfile)

            if any(strcmp(stadat.chans.component,'Z'))
                eqar(is).datZ = interp1(stadat.tt,stadat.dat(:,strcmp(data.chans.component,'Z')),eqar(is).tt);
            end
            if any(strcmp(data.chans.component,'N')) & any(strcmp(data.chans.component,'E'))
                datN = interp1(stadat.tt,stadat.dat(:,strcmp(stadat.chans.component,'N')),eqar(is).tt);
                datE = interp1(stadat.tt,stadat.dat(:,strcmp(stadat.chans.component,'E')),eqar(is).tt);

                foraz = mod(stadat.seaz+180,360);
                eqar(is).datR =  datN*sind(foraz) + datE*sind(foraz);
                eqar(is).datT = -datN*sind(foraz) + datE*cosd(foraz);
            end
            
            if isfield(datinfo,'xcor')
            if datinfo(is).xcor
                eqar(is).dT = [];
                eqar(is).acor = [];
                eqar(is).abs_arrT = [];
                eqar(is).par_dT = [];
            end
            end
            if isfield(datinfo,'dtstar')
            if datinfo(is).dtstar
                eqar(is).snr_wf = [];
                eqar(is).frq = [];
                eqar(is).specn = [];
                eqar(is).specs = [];
                eqar(is).specss = [];
                eqar(is).fcross = [];
                eqar(is).dtstar = [];
                eqar(is).std_dtstar = [];
                eqar(is).par_dtstar_specR = [];
            end
            end
            
            if ifsave
                save(eqarfile,'eqar');
            end

        end % loop on eqarchans
        
        datinfo(is).NEZ    = true;
        datinfo(is).rmresp = true;
        datinfo(is).rmtilt = false;
        datinfo(is).rmcomp = false;
        datinfo(is).xcor = false;
        datinfo(is).dtstar = false;
        fprintf(' > ')
        if ifsave
            save(phdatinfofile,'datinfo');
        end
        
    end % loop on phases
    
    fprintf('fixed\n')
end % loop on events
