% for each event, cycle through data and build structures with all
% necessary data for each arrival of interest. 
clear all
% close all
addpath('matguts')

%% parameters
overwrite = true;
phase = 'S';
resamprate = 5 ; % new, common sample rate
wind = [-200 200]; % seconds before and after arrival to save data for this arrival


%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% DATA DIRECTORY (top level)
datadir = '/Volumes/DATA_mini2/CASCADIA/DATA/'; % needs final slash

%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

% GET EVENTS DATA
db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,elats,elons,edeps,evtimes,mags] = dbgetv(dbor,'orid','lat','lon','depth','time','ms');
norids = dbnrecs(dbor);
dbclose(db);

for ie = 1:270 % 44:norids % loop on orids
    fprintf('\n Orid %.0f %s \n\n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
    evdir = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
    if ~exist([datadir,evdir,'_datinfo.mat'],'file'), fprintf('No data at all for this event\n'), continue, end
    datinfofile = [datadir,evdir,'_datinfo'];
    load(datinfofile)
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end
    nstas = length(datinfo);

    ofile = [datadir,evdir,'_EQAR_',phase];
    if exist([ofile,'_Z.mat'],'file')
        if overwrite
            yn = 'y';
        else
            yn = input(sprintf('%s EQAR file exists, overwrite? [y/n] ',phase),'s');
        end
        if strcmp(yn,'y')
            delete([ofile,'_Z.mat'])
            delete([ofile,'_R.mat'])
            delete([ofile,'_T.mat'])
        else
            fprintf('ok, skipping\n')
            continue
        end
    end
    
    % RESULTS STRUCTURE
    eqar = struct('phase',phase,...
                  'sta',{datinfo.sta}','slat',[],'slon',[],'selev',[],...
                  'gcarc',0,'seaz',0,'rayp',0,'pred_arrT',0,...
                  'tt',[],'datZ',[],'datR',[],'datT',[],'datH',[],'corZ',[],'samprate',resamprate);

    wlen = diff(wind)*resamprate; % window length in samples
    
    % LOOP ON STAS
    for is = 1:nstas
        fprintf('Station %-8s ',sprintf('%s...',datinfo(is).sta))
        % APPLY DATA QUALITY CONDITIONS
%         if datinfo(is).crap == 1; 
%             fprintf('crap, skipping\n'),continue 
%         end
        if ~datinfo(is).rmresp
            fprintf('must have response removed\n'), continue
        end
        
        % GET STATION + ARRIVAL INFO
        load([datadir,evdir,datinfo(is).sta,'.mat']); % load sta data for this evt
        
        if ~any(strcmp({data.phases.phase},phase))
            fprintf('No %s arrival at this sta\n',phase), continue
        end
        
        % station details
        statmp = struct2cell(data.station);
        [eqar(is).sta,eqar(is).slat,eqar(is).slon,eqar(is).selev] = deal(statmp{:});
        
        % station-event details
        eqar(is).gcarc = data.gcarc;
        eqar(is).seaz = data.seaz;

        % arrival details
        ip = find(strcmp({data.phases.phase},phase),1,'first');
        eqar(is).rayp = data.phases(ip).rayparameter;
        eqar(is).pred_arrT = data.phases(ip).artime;
        
        % check resampling
        if resamprate>data.samprate
            error('resamprate would alias data - use smaller resamprate')
        end
            
        % GET DATA
        tt = [data.phases(ip).artime + wind(1):1./resamprate:data.phases(ip).artime + wind(2)];
        tt = tt(1:wlen);
        eqar(is).tt = tt;
        if any(strcmp(data.chans.component,'Z'))
            eqar(is).datZ = interp1(data.tt,data.dat(:,strcmp(data.chans.component,'Z')),tt);
        end
        if any(strcmp(data.chans.component,'H'))
            eqar(is).datH = interp1(data.tt,data.dat(:,strcmp(data.chans.component,'H')),tt);
        end
        
        if ~datinfo(is).NEZ
            fprintf(' not rotated\n'), continue
        else
            if data.chans.azimuth(strcmp(data.chans.component,'N')) ~=0, error('Bad instrument orientation\n'), end
            datN = interp1(data.tt,data.dat(:,strcmp(data.chans.component,'N')),tt);
            datE = interp1(data.tt,data.dat(:,strcmp(data.chans.component,'E')),tt);
            
            foraz = mod(data.seaz+180,360);
            eqar(is).datR =  datN*sind(foraz) + datE*sind(foraz);
            eqar(is).datT = -datN*sind(foraz) + datE*cosd(foraz);
        end
        
        % GET NOISE DATA IF IT EXISTS
        if data.rmtilt || data.rmcomp
            if data.noise_cor.tt(1)==-100, data.noise_cor.tt = data.noise_cor.tt+evtimes(ie); end % correct if in event rather than absolute time ref frame
            eqar(is).corZ = interp1(data.noise_cor.tt,data.noise_cor.Zdat,tt);
            if any(isnan(eqar(is).corZ)), error('Bad interp?!'); end
        end
        
        fprintf('got data\n')    
    end % loop on stas

    % SAVE
    save([ofile,'_Z.mat'],'eqar')
    save([ofile,'_R.mat'],'eqar')
    save([ofile,'_T.mat'],'eqar')
    
    copyfile([datinfofile,'.mat'],[datinfofile,'_',phase,'.mat'])

    
end% loop on orids
