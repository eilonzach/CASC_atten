function [ Qs,Qp,Vs,Vp ] = QV_at_z( age,Z,gs,frq, rho )
% [ Qs,Qp ] = Qs_z( age )
%   Function to calculate a Q(z) profile, assuming a certain geotherm and
%   then using the Jackson and Faul, 2010 treatment to predict Qs etc.

if nargin < 2
    Z = 10:5:200;
end
if nargin < 3
    gs = 0.001;
end
if nargin < 4
    frq = 1;
end
if nargin < 5
    rho = 3300;
end

TempZ = geotherm( age,'plate',Z,1350);

addpath('/Users/zeilon/Documents/MATLAB/geoff/VpVsQ')

for ii = 1:length(Z)
[qinv(ii),gg(ii),ks(ii)]=fjcalc(TempZ(ii), gs, frq, Z(ii)/30);
end

Qs = 1./qinv;
Qp = (9/4)*Qs; % using classic relationship
Vs = sqrt(gg*1e9/rho);
Vp = sqrt((ks + 1.333*gg)*1e9/rho);

% plot(Qs,-Z)
% xlim([0 500])



end

