% script to go through data and rotate each station's channels 
% system from whatever the OBS values were to the N,E,Z
clear all

overwrite = true;

% path to orientation data
oriedir = '~/Work/CASCADIA/ORIENTATIONS_RF/';

% antelope db details
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascattendb';

% path to top level of directory tree for data
datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash

%% get to work

oriefiles = dir(oriedir); 

db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,elats,elons,edeps,evtimes,mags] = dbgetv(dbor,'orid','lat','lon','depth','time','ms');
norids = dbnrecs(dbor);
dbclose(db);

for ie = 85:85 % 1:norids % loop on orids
    fprintf('\n Orid %.0f %s \n\n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
    evdir = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
%     stafiles = dir([datadir,evdir]); stafiles = stafiles(3:end);
    load([datadir,evdir,'_datinfo.mat'])
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end

    for is = 1:length(datinfo) % loop on stas
        if datinfo(is).NEZ, continue, end % skip if already rotated
        
        sta = datinfo(is).sta; % sta name
        fprintf('Station %s...',sta)
        load([datadir,evdir,sta,'.mat']); % load sta data for this evt
        if strcmp(sta,'M02CO') || strcmp(sta,'M04CO'), sta = sta(1:4); end % rename problem stations
        
        yesE = any(strcmp(datinfo(is).chans,'E'));
        yesN = any(strcmp(datinfo(is).chans,'N'));
        yes1 = any(strcmp(datinfo(is).chans,'1'));
        yes2 = any(strcmp(datinfo(is).chans,'2'));
        
         % CASE: NO HORIZONTALS
        if (~yesE && ~yesN && ~yes1 && ~yes2)
            fprintf(' have no horizontals!\n')
        end
        
        % CASE: ALREADY NEZ
        if (yesE && yesN)
            fprintf(' no need to rotate, making sure NEZ order.\n')
            chans = data.chans;
            ich = find(strcmp(chans.component,'H'));
            ice = find(strcmp(chans.component,'E'));
            icn = find(strcmp(chans.component,'N'));
            icz = find(strcmp(chans.component,'Z'));
            
            newdat = data.dat(:,[ich,icn,ice,icz]);
            chans.name = chans.name([ich,icn,ice,icz]);
            chans.component = chans.component([ich,icn,ice,icz]);
            chans.azimuth = chans.azimuth([ich,icn,ice,icz]);
             % put back into data structure
            data.dat   = newdat;
            data.chans = chans;
            data.NEZ = true;
            
            datinfo(is).NEZ = true;                            %#ok<SAGROW>
            datinfo(is).chans = chans.component;               %#ok<SAGROW>
        end
        
        % CASE: Have E or N but neither 1 or 2
        if ((yesE && ~yesN) || (~yesE && yesN)) && ~yes1 && ~yes2
            fprintf(' have one N or E horiz but not the other; no need to rotate.\n')
            datinfo(is).NEZ = true;                            %#ok<SAGROW>
            data.NEZ = true;
        end
        
        % CASE: Have 1 or 2 but neither E or N
        if ((yes1 && ~yes2) || (~yes1 && yes2)) && ~yesN && ~yesE
            fprintf(' have one 1 or 2 horiz but not the other - cannot rotate\n')
        end
        
        % CASE: CHANS 1&2 ==> ROTATE
        if (yes1 && yes2)
            fprintf(' rotating 1,2 ==> N,E...')
            
            chans = data.chans;
            ich = find(strcmp(chans.component,'H'));
            ic1 = find(strcmp(chans.component,'1')); % instrument N
            ic2 = find(strcmp(chans.component,'2')); % instrument E
            icz = find(strcmp(chans.component,'Z'));
            
            if any(chans.azimuth([ic1,ic2]))
                orientation = chans.azimuth(ic1);
                fprintf(' using IRIS orientation')
            elseif ~any(chans.azimuth([ic1,ic2]))
                orfile = dir([oriedir,sta,'*']);
                if isempty(orfile), fprintf(' no orientation given, CANNOT ROTATE\n'), continue, end
                or = load([oriedir,orfile.name]);
                orientation = or.orientation;
                fprintf(' using Helen RF orientation')
            end

            newdat = zeros(size(data.dat));
            sor = sind(orientation);
            cor = cosd(orientation);

            newdat(:,ic1) = cor*data.dat(:,ic1) - sor*data.dat(:,ic2);    % NORTH
            newdat(:,ic2) = sor*data.dat(:,ic1) + cor*data.dat(:,ic2);    % EAST
            newdat(:,ich) = data.dat(:,ich);
            newdat(:,icz) = data.dat(:,icz);
            
            % alter chans
            chans.name([ic1,ic2]) = {[chans.name{ic1}(1:2),'N'],[chans.name{ic2}(1:2),'E']};
            chans.component([ic1,ic2]) = {'N','E'};
            chans.azimuth([ic1,ic2]) = orientation +[0 90];
    
            % put back into data structure
            data.dat   = newdat;
            data.chans = chans;
            data.NEZ = true;
            
            % save            
            fprintf(' - done rotating.\n')
            datinfo(is).NEZ = true;                            %#ok<SAGROW>
            datinfo(is).chans = chans.component;               %#ok<SAGROW>
        end
        
        save([datadir,evdir,sta],'data')

    end % loop on stas
    save([datadir,evdir,'_datinfo'],'datinfo')
    
	%% sum up
    fprintf(' STA  CHAN  NEZ  resp  tilt  comp\n')
    for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f     %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).rmtilt,datinfo(is).rmcomp); end

end
