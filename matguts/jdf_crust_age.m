function [ ages,chrons,F ] = jdf_crust_age( lats,lons )
% [ ages,chrons,F ] = jdf_crust_age( lats,lons )
%  Function to take lat,lon points and interpolate a grid of seafloor ages
%  to compute the age of the crust at the points specified, on the Juan de
%  Fuca or Gorda plates

%% key for chron name ===> age
c2afile = '~/Work/CASCADIA/DATA/mapdata/isochrons_dwilson/JGR1993/chron_ages_ze.txt';
fid = fopen(c2afile,'r');
A = textscan(fid,'%s %f','CommentStyle','#');
fclose(fid);
chron.name = A{1};
chron.age = A{2};

%% read isochron file
chronlnlt_file = '~/Work/CASCADIA/DATA/mapdata/isochrons_dwilson/JGR1993/isochrons.lnlt';

chlat = [];
chlon = [];
chage = [];

fid = fopen(chronlnlt_file,'r');
eof = 0;
count = 1;
while eof==0
    str = fgets(fid);
    if strcmp(str(1),'>') % read in segment separators
        
        B = textscan(str,'> %s %s');
        
        if strcmp(B{2},''),     % segment separator
            chlat(count,1) = nan;
            chlon(count,1) = nan;
            chage(count,1) = nan;
            count = count+1;
        else                    % New age group
            cnm = B{2};
            cnm = strtok(cnm,'_j');
            if ~any(strcmp(chron.name,cnm)), error('this chron has no defined age!'), end
            age = chron.age(strcmp(chron.name,cnm));   
        end
    else % lat/lon datum
        C = textscan(str,'%f %f');
        chlon(count,1) = C{1};
        chlat(count,1) = C{2};
        chage(count,1) = age;
        count = count+1;
    end

    eof = feof(fid);
end
fclose(fid);

%% read ridge axis file
rdglnlt_file = '~/Work/CASCADIA/DATA/mapdata/isochrons_dwilson/JGR1993/ridge.lnlt';

fid = fopen(rdglnlt_file,'r');
eof = 0;
while eof==0
    str = fgets(fid);
    if strcmp(str(1),'>') % segment separator  
        chlat(count,1) = nan;
        chlon(count,1) = nan;
        chage(count,1) = nan;
        count = count+1;
    else % lat/lon datum
        C = textscan(str,'%f %f');
        chlon(count,1) = C{1};
        chlat(count,1) = C{2};
        chage(count,1) = 0;
        count = count+1;
    end
    eof = feof(fid);
end
fclose(fid);




% figure(1),clf, hold on
% scatter(chlon(chage<15),chlat(chage<15),60,chage(chage<15), 'filled')

nnan = ~isnan(chage);
F = scatteredInterpolant(chlon(nnan),chlat(nnan),chage(nnan),'linear','none');

ages = F(lons,lats);
ages(ages > 1.5*max(chage(nnan))) = nan;

chrons = struct('lat',chlat,'lon',chlon,'age',chage);

end

