function [ dist_deg,dist_km,jdf_or_gor ] = dist_to_rdg( lats,lons )
% [ dist_km,dist_deg ] = dist_to_rdg( lats,lons )
%   Quick and dirty function to find the distance (East) from the JdF or
%   Gorda ridges - just using approximate lines for the ridge axis and
%   working out distance to these lines

lats = lats(:); lons = lons(:);

dist_deg = zeros(size(lats));
dist_km = zeros(size(lats));

kpd = 78; % km per degree of longitude

%% approximate lines of the ridges
roughjdf = [-130.5,44.36;-128.7,48.98];
roughgor = [-127.6,40.49;-126.6,43.01];

roughblanco = [-126.6,43.01;-130.4,44.36];% in order to parse if sta is on JdF or Gor plate

% work out distance to each ridge/transform
Xrdg_jdf = dist2line(roughjdf(1,:),roughjdf(2,:),[lons,lats]); % positive means E of ridge
Xrdg_gor = dist2line(roughgor(1,:),roughgor(2,:),[lons,lats]);  % positive means E of ridge
Xblanco = dist2line(roughblanco(1,:),roughblanco(2,:),[lons,lats]);  % positive means N of transform

jdf_or_gor = Xblanco > 0; 

dist_deg(jdf_or_gor) = Xrdg_jdf(jdf_or_gor);
dist_deg(~jdf_or_gor) = Xrdg_gor(~jdf_or_gor);

dist_km = kpd*dist_deg;

end

