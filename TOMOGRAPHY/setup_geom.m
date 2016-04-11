function [ data,par ] = setup_geom( data,par )
% [ data,par ] = setup_geom( data,par )
% setup the geometry of the model space for the inverse problem.



%% ========== define model space ==========================================
[stax,stay] = project_xy( par, data.stn.lat, data.stn.lon );
data.stn.stax = stax;
data.stn.stay = stay;

% increase node spacing from centre out and down
maxx = round_level(max(stax),par.dh); 
minx = round_level(min(stax),par.dh);
maxy = round_level(max(stay),par.dh); 
miny = round_level(min(stay),par.dh);
% lx = maxx-minx;
% ly = maxy-miny;
xx = [minx:par.dh:maxx];
yy = [miny:par.dh:maxy];

% endy = cumsum(round_level(linspace(par.dh,par.dh_max,ceil(lx/(mean([par.dh,par.dh_max])))),5));
endxy = cumsum(round_level(linspace(par.dh,par.dh_max,ceil(par.zmax/(mean([par.dh,par.dh_max])))),5));

% define nodes
xx = [minx-fliplr(endxy) xx maxx+endxy]';
yy = [miny-fliplr(endxy) yy maxy+endxy]';
zz = [0 par.zmin par.zmin+cumsum(round_level(linspace(par.dz,par.dz_max,floor((par.zmax-par.zmin)/(mean([par.dz,par.dz_max])))),5))]';
% put nodes into vectors
[X,Y,Z] = meshgrid(xx,yy,zz);
    % cycle over y then x then z.
mx = reshape(X,numel(X),1);
my = reshape(Y,numel(X),1);
mz = reshape(Z,numel(X),1);
[mlt, mln] = project_xy(par, mx, my, 'inverse');

dxx = 0.5*([0;diff(xx)]+[diff(xx);0]);
dyy = 0.5*([0;diff(yy)]+[diff(yy);0]);
% dzz = 0.5*([0;diff(zz)]+[diff(zz);0])
dzz = [0;diff(zz)];
[X,Y,Z] = meshgrid(dxx,dyy,dzz);
mdx = reshape(X,numel(X),1);
mdy = reshape(Y,numel(X),1);
mdz = reshape(Z,numel(X),1);

%% save into par structure
par.mx = mx;
par.my = my;
par.mz = mz;
par.mlt = mlt;
par.mln = mln;
par.mdx = mdx;
par.mdy = mdy;
par.mdz = mdz;

par.xx = xx;
par.yy = yy;
par.zz = zz;

par.nx = length(xx);
par.ny = length(yy);
par.nz = length(zz);
par.nmodel = par.nx*par.ny*par.nz;

par.min_x = min(xx);
par.max_x = max(xx);
par.min_y = min(yy);
par.max_y = max(yy);
par.min_z = min(zz);
par.max_z = max(zz);

%% average velocities 
% calculated as mean over half-bins above and below par.zz
par.vz = zeros(size(par.zz));
zz = [0;par.zz;par.zz(end)];
for iz = 1:length(par.vz)
    zbo = 0.5*(zz(iz+1)+zz(iz+2));
    zto = 0.5*(zz(iz)+zz(iz+1));
    par.vz(iz) = quad(@(z)vel_profile(par.PS,z),zto,zbo)/(zbo-zto);    
end


par.mvav   = interp1(par.zz,par.vz,par.mz);
% [~,par] = make_start_model(par); % 

end
% 
% function V = vel_profile_integrand(z)
%     V = vel_profile(par.PS,z);
% end