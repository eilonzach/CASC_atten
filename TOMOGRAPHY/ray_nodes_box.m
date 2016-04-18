function [n_indx, n_lens ] = ray_nodes_box(rx,ry,rz,dr,par)
%ray x,y,z,dr  in km

Nr = length(dr);

n_indx = nan(8*Nr,1);
n_lens = nan(8*Nr,1);

xind = ([1:par.nx]-1)*par.ny;
yind =  [1:par.ny];
zind = ([1:par.nz]-1)*par.nx*par.ny;

for ir = 1:Nr
    xlo = max( par.xx( par.xx <= rx(ir) ) );
    xhi = min( par.xx( par.xx >  rx(ir) ) );
    ylo = max( par.yy( par.yy <= ry(ir) ) );
    yhi = min( par.yy( par.yy >  ry(ir) ) );
    zmi = max( par.zz( par.zz <= rz(ir) ) ); 
        if isempty(zmi), zmi = min(par.zz); end % for points above surf
    zma = min( par.zz( par.zz >  rz(ir) ) );
        if isempty(zma), zma = max(par.zz); end % for points below base
    
    indx = zeros(8,1);
    try
    indx(1) = xind(par.xx==xlo) + yind(par.yy==ylo) + zind(par.zz==zmi) ; % 000
    indx(2) = xind(par.xx==xhi) + yind(par.yy==ylo) + zind(par.zz==zmi) ; % 100
    indx(3) = xind(par.xx==xlo) + yind(par.yy==yhi) + zind(par.zz==zmi) ; % 010
    indx(4) = xind(par.xx==xlo) + yind(par.yy==ylo) + zind(par.zz==zma) ; % 001
    indx(5) = xind(par.xx==xhi) + yind(par.yy==yhi) + zind(par.zz==zmi) ; % 110
    indx(6) = xind(par.xx==xhi) + yind(par.yy==ylo) + zind(par.zz==zma) ; % 101
    indx(7) = xind(par.xx==xlo) + yind(par.yy==yhi) + zind(par.zz==zma) ; % 011
    indx(8) = xind(par.xx==xhi) + yind(par.yy==yhi) + zind(par.zz==zma) ; % 111
    catch
        keyboard
    end
    
    fx = 1 - (rx(ir)-xlo)/(xhi-xlo); % fraction between xlo and xhi (0-1)
    fy = 1 - (ry(ir)-ylo)/(yhi-ylo); % fraction between ylo and yhi (0-1)
    fz = 1 - (rz(ir)-zmi)/(zma-zmi); % fraction between zmi and zma (0-1)
    
    if zma==zmi, fz = 0; end % in case node at surface, avoid Inf
    
    vals = [fx*fy*fz;
            (1-fx)*fy*fz;
            fx*(1-fy)*fz;
            fx*fy*(1-fz);
            (1-fx)*(1-fy)*fz;
            (1-fx)*fy*(1-fz);
            fx*(1-fy)*(1-fz);
            (1-fx)*(1-fy)*(1-fz)]; % linearly interpolate in 3D
     
%     if abs(sum(vals)-1) > 1e-15 % check that vals sum to 1
%         error('Summed sensitivities to box corners should be 1'); 
%     end
     
    n_indx([1:8] + 8*(ir-1)) = indx;
    n_lens([1:8] + 8*(ir-1)) = vals .* dr(ir);
    
end

if any(isnan(n_indx)), error('Some node indices not assigned'), end
if any(isnan(n_lens)), error('Some node values not assigned'), end

% now have to group lengths for each node - as n_indx will contain several
% repetitions
un_indx = unique(n_indx); % find unique node indices
un_lens = zeros(size(un_indx));
for n = 1:length(un_indx) % loop over unique nodes
    un_lens(n) = sum(n_lens(n_indx==un_indx(n))); % summed sensitivity is sum of lengths thru that node
end

% rename
n_indx = un_indx;
n_lens = un_lens;
    
end

