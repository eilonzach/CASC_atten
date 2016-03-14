function [rho,KS,G,VSh,VSv,VPh,VPv,Z ] = HeFESTo_getval( Ti, Pi )
% [rho,KS,G,VSh,VSv,VPh,VPv,Z ] = HeFESTo_getval( Ti, Pi )
% 
% INPUT
%  Ti, Pi are vectors of temp. and pressure of desired points, in C / GPa

ifile = 'HeFESTo.out';
[P,Z,T,rho,KS,G,VBh,VSh,VPh,VBr,VSr,VPr,VBv,VSv,VPv] = HeFESTo_read(ifile);

% find appropriate row for each of Ti,Pi
M = length(P);
N = length(Ti);
P = P*ones(1,N); Pi = ones(M,1)*Pi(:)';
T = T*ones(1,N); Ti = ones(M,1)*Ti(:)';
dP = abs(P - Pi);
dT = abs(T - Ti);
dE = dP./mean(mean(P)) + dT./mean(mean(T));

[~,indx] = min(dE);
indx = indx(:);

Z   = Z(indx);
rho = rho(indx)*1000;
KS  = KS(indx);
G   = G(indx);
VSh = VSh(indx);
VPh = VPh(indx);
VSv = VSv(indx);
VPv = VPv(indx);


end

