function [rx,ry,rz,Dr] = ray_trace_1D(PS, p, baz, elat,elon,slat, slon,stn_h, par)

% determine if potentially turning in model
max_depth = par.max_z;
vdz = 5;
rldr = par.kidd;

% if not turn, else consider turn
zz1 = [0:vdz:(max_depth-vdz)]';
zz2 = [vdz:vdz:max_depth]';
v1  = vel_profile(PS,zz1,par.Vavmod,stn_h);
v2  = vel_profile(PS,zz2,par.Vavmod,stn_h);
r1 = 6371-zz1;
r2 = 6371-zz2;
zf1 = 6371*log(6371./r1);
zf2 = 6371*log(6371./r2);
vf1 = (6371./r1).*v1;
vf2 = (6371./r2).*v2;
u1 = 1./vf1;
u2 = 1./vf2; 
dv = vf2-vf1;
dz = zf2-zf1;
b = dv ./ dz;
const_indx = (b == 0);

X = zeros( length(v1), 1 );
X(const_indx) = constv_dist( u1(const_indx), dz(const_indx), p );
X(~const_indx) = gradv_dist(b(~const_indx),u1(~const_indx),u2(~const_indx),p);
X = [0; X];
X(imag(X)>0) = NaN;
X = X(~isnan(X));
Xd = cumsum(X);
Xd = r2d(Xd/6371); % convert km to degrees
[rlat,rlon] = reckon( slat, slon, Xd, baz );
[rayx,rayy] = project_xy( par, rlat, rlon );
rayz = [zz1(1); zz2]; 
rayz = rayz(1:length(X));
dz = dz(1:length(X)-1);

% calc Dr
r_scale = (1/6371).*(6371-rayz);
X2 = r_scale.*X;
Dr1 = (X2(2:end).^2 + dz.^2).^(1/2);
Dr1 = [0; Dr1];
Dr2 = cumsum(Dr1);
Dr = [0:rldr:max(Dr2)]';
Dr = [Dr; max(Dr2)];
% cDr = Dr;
rx = interp1(Dr2,rayx,Dr);
ry = interp1(Dr2,rayy,Dr);
rz = interp1(Dr2,rayz,Dr);
% nn = 1:round(length(Dr)/length(rz)):length(Dr)+10;
% nn = nn(1:length(rz));
% Dr = interp1(1:length(Dr),Dr,nn)'; % this will be exactly the same Dr as above
% cDr = interp1(1:length(cDr),cDr,nn)'; % this is the same as Dr

if min(rx)<par.min_x
    ind_end = find(rx<(par.min_x),1,'first');
    rx = rx(1:(ind_end));
    ry = ry(1:(ind_end));
    rz = rz(1:(ind_end));
    Dr = Dr(1:(ind_end));
%     cDr = cDr(1:(ind_end));
end
if max(rx)>par.max_x
    ind_end = find(rx>(par.max_x),1,'first');
    rx = rx(1:(ind_end));
    ry = ry(1:(ind_end));
    rz = rz(1:(ind_end));
    Dr = Dr(1:(ind_end));
%     cDr = cDr(1:(ind_end));
end
if min(ry)<par.min_y
    ind_end = find(ry<(par.min_y),1,'first');
    rx = rx(1:(ind_end));
    ry = ry(1:(ind_end));
    rz = rz(1:(ind_end));
    Dr = Dr(1:(ind_end));
%     cDr = cDr(1:(ind_end));
end
if max(ry)>par.max_y
    ind_end = find(ry>(par.max_y),1,'first');
    rx = rx(1:(ind_end));
    ry = ry(1:(ind_end));
    rz = rz(1:(ind_end));
    Dr = Dr(1:(ind_end));
%     cDr = cDr(1:(ind_end));
end


%-------------------- sub-functions --------------------------------------
%
function x = constv_dist( u, dz, p )

x = (p * dz) ./ eta(u,p);


function x = gradv_dist( b, u1, u2, p )

x1 = eta(u1,p) ./ (b.*u1*p);
x2 = eta(u2,p) ./ (b.*u2*p);

x = x1 - x2;

%
function val = eta( u, p )

val = sqrt( u.^2 - p.^2 );
