function [n_indx, n_vals ] = ray_nodes_ze(rx,ry,rz,dr,rmax,par,cf,gcarc,time)
%rayx,y,z,dr  in km
% first pass - find nodes nearby ray
l_r = length(dr);
kdr = dr(2)-dr(1); % path length element
nn = 2:4:l_r-1;
rayxyz2 = [ rx(nn) ry(nn) rz(nn) ];
rmax2 = rmax(nn);
 
nn=1;
n_indx = zeros(20000,1);
for qq = 1:length(rayxyz2(:,1))
    dists = ( (par.mx - rayxyz2(qq,1)).^2 ...
            + (par.my - rayxyz2(qq,2)).^2 ...
            + (par.mz - rayxyz2(qq,3)).^2 ).^(1/2);
%     rmax_eff = (rayxyz2(qq,3)/par.max_z).^0.5 * rmax;    
    rmax_z = sqrt( (4*kdr)^2 + rmax2(qq)^2); % for a first pass using every fourth node
    ind1 = find(dists<rmax_z);
    n_indx(nn:(nn+length(ind1)-1)) = ind1;
    nn = nn+length(ind1);
end

n_indx = n_indx(n_indx>0);
n_indx = unique(n_indx); % these indices are the numbers of the nodes within the sensitivity 

l_nindx = length(n_indx);

% calc Rn for all nearby nodes
r_scale = (6371-rz)/6371;

r_scale = r_scale*ones(1,l_nindx);
rayx = rx*ones(1,l_nindx);
rayy = ry*ones(1,l_nindx);
rayz = rz*ones(1,l_nindx);
nodx = ones(l_r,1)*par.mx(n_indx)';
nody = ones(l_r,1)*par.my(n_indx)';
nodz = ones(l_r,1)*par.mz(n_indx)';

dists = ( (r_scale.*(rayx-nodx)).^2 ...
        + (r_scale.*(rayy-nody)).^2 ...
        + (          rayz-nodz ).^2 ).^(1/2);
[n_rn,ind1] = min(dists); 
n_rn = n_rn';
n_rf = rmax(ind1);

% n_rn2 = zeros(l_nindx,1);
% n_rf2 = zeros(l_nindx,1);
% 
% for qq = 1:length(n_indx)
%    dists = ( (r_scale.*(rx-par.mx(n_indx(qq)))).^2 ...
%            + (r_scale.*(ry-par.my(n_indx(qq)))).^2 ...
%            + (rz -par.mz(n_indx(qq))).^2 ).^(1/2);   
% 	[n_rn2(qq),ind1] = min(dists);
% 	n_rf2(qq) = rmax(ind1);
% end

ind2 = n_rn>n_rf;
n_indx(ind2)= [];
n_rn(ind2)  = [];
n_rf(ind2)  = [];

ind3 = n_indx <= par.nx*par.ny;
n_indx(ind3)= [];
n_rn(ind3)  = [];
n_rf(ind3)  = [];

% vols = par.mdx(n_indx).*par.mdy(n_indx).*par.mdz(n_indx);
vavs = par.mvav(n_indx);
% function to calc sensitivity vals
[ n_indx, n_vals ] = f1_ker(n_indx,n_rn,n_rf,vavs,time);

end

