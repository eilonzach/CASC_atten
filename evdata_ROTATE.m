% script to go through data and rotate each station's channels 
% system from whatever the OBS values were to the N,E,Z
clear all

overwrite = true;

% path to orientation data
oriedir = '~/Work/CASCADIA/ORIENTATIONS_RF/';

% antelope db details
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';

% path to top level of directory tree for data
datadir = '/VOLUMES/DATA_mini/CASCADIA/DATA/'; % needs final slash

%% get to work

oriefiles = dir(oriedir); 

db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,elats,elons,edeps,evtimes,mags] = dbgetv(dbor,'orid','lat','lon','depth','time','ms');
norids = dbnrecs(dbor);
dbclose(db);

for ie = 215:215 % 1:norids % loop on orids
    fprintf('\n Orid %.0f %s \n\n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
    evdir = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
    datinfofile = [datadir,evdir,'_datinfo'];

    if ~exist([datinfofile,'.mat'],'file'), fprintf('No data at all for this event\n'), continue, end
    load(datinfofile)
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end

    for is = 1:length(datinfo) % loop on stas
        if datinfo(is).NEZ, continue, end % skip if already rotated
        
        sta = datinfo(is).sta; % sta name
        fprintf('Station %s...',sta)
        load([datadir,evdir,sta,'.mat']); % load sta data for this evt
        
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
            fprintf('maybe minor rotation + making sure NEZ order.\n')
            chans = data.chans;
            ich = find(strcmp(chans.component,'H'));
            ice = find(strcmp(chans.component,'E'));
            icn = find(strcmp(chans.component,'N'));
            icz = find(strcmp(chans.component,'Z'));
            
            newdat = zeros(size(data.dat));
            sor = sind(chans.azimuth(icn));
            cor = cosd(chans.azimuth(icn));
            
            ic_new = [ich icn ice icz];

            newdat(:,ich) = data.dat(:,ich);
            newdat(:,icn) = cor*data.dat(:,icn) - sor*data.dat(:,ice);    % NORTH
            newdat(:,ice) = sor*data.dat(:,icn) + cor*data.dat(:,ice);    % EAST
            newdat(:,icz) = data.dat(:,icz);
            newdat = newdat(:,ic_new);
            
            chans.name = chans.name([ich,icn,ice,icz]);
            chans.component = chans.component([ich,icn,ice,icz]);
            chans.azimuth = [chans.azimuth(ich), 0 90 chans.azimuth(icz)];
            
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
            
            ic_new = [ich ic1 ic2 icz];

            orientation = [];
            if any(chans.azimuth([ic1,ic2]))
                orientation = chans.azimuth(ic1);
                if chans.azimuth(ic2)-chans.azimuth(ic1) ~= 90, fprintf('Channels in unexpected order! Skipping...\n'), continue, end
                fprintf(' using IRIS orientation')
            elseif ~any(chans.azimuth([ic1,ic2]))        
%                 if strcmp(sta,'M02CO') || strcmp(sta,'M04CO'), sta_use = sta(1:4); else sta_use = sta; end % rename problem stations
%                 orfile = dir([oriedir,sta_use,'*']);
%                 if ~isempty(orfile) % use helen's orientation
%                     or = load([oriedir,orfile.name]);
%                     orientation = or.orientation;
%                     fprintf(' using Helen RF orientation')
%                 else
                db = dbopen([dbdir,dbnam],'r');
                dbsch = dblookup_table(db,'sitechan');
                dbsch.record = dbfind(dbsch,sprintf('sta == "%s" && chan == "%s"',sta,chans.name{ic1}));
                [orientation1] = dbgetv(dbsch,'hang'); % hang of nominal north
                dbsch.record = dbfind(dbsch,sprintf('sta == "%s" && chan == "%s"',sta,chans.name{ic2}));
                [orientation2] = dbgetv(dbsch,'hang'); % hang of nominal east                
                dbclose(db);
                if mod(orientation1+90,360) == orientation2 % check north is 90 acw of east
                    orientation = orientation1;
                    fprintf(' using OBSIP orientation')
                end
%                 end
            end
            if isempty(orientation)
                fprintf(' no orientation given, CANNOT ROTATE\n'), continue
            end
            
                
            newdat = zeros(size(data.dat));
            sor = sind(orientation);
            cor = cosd(orientation);

            newdat(:,ic1) = cor*data.dat(:,ic1) - sor*data.dat(:,ic2);    % NORTH
            newdat(:,ic2) = sor*data.dat(:,ic1) + cor*data.dat(:,ic2);    % EAST
            newdat(:,ich) = data.dat(:,ich);
            newdat(:,icz) = data.dat(:,icz);
            
            newdat = newdat(:,ic_new);

            % alter chans
            chans.name([ic1,ic2]) = {[chans.name{ic1}(1:2),'N'],[chans.name{ic2}(1:2),'E']};
            chans.component([ic1,ic2]) = {'N','E'};
            chans.azimuth([ic1,ic2]) = [0 90];
    
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
    save(datinfofile,'datinfo')

    copyfile([datinfofile,'.mat'],[datinfofile,'_P.mat'])
    copyfile([datinfofile,'.mat'],[datinfofile,'_S.mat'])
    
	%% sum up
    fprintf(' STA  CHAN  NEZ  resp  tilt  comp\n')
    for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f     %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).rmtilt,datinfo(is).rmcomp); end

end
