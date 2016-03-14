function [ dtstar,dT,A0,misfit, E ] = invert_Aphi_4_dtdtstar( As,phis,freqs, wts,amp2phiwt)
%[ dtstar,dT,A0 ] = invert_Aphi_4_dtdtstar( As,phis,freqs, wts)
%   Detailed explanation goes here
% A = A0*exp(-pi*dtstar*f)
% ln(A) = ln(A0) - pi*dtstar*f 
% phi = (ln(f) - ln(fNq))*dtstar/pi + dT

if nargin < 4
    wts = ones(size(As));
end
if nargin < 5
    amp2phiwt = 1;
end
As = As(:); phis = phis(:); wts = wts(:);

N = length(As);

% Model parameter vector will be [ln(A0); dtstar; dT]

% data vector - amplitudes and phis
d = [log(As);phis];
% G matrix - relationship between 
G = zeros(2*N,3);
G([1:N],1) = 1;
G([1:N],2) = -pi*freqs;
G(N+[1:N],2) = 1 - log(freqs)./pi;
G(N+[1:N],3) = 1;

% weight
wts = [amp2phiwt*wts;wts];
dw = diag(wts)*d;
Gw = diag(wts)*G;

m = (Gw'*Gw)\Gw'*dw;    

dtstar = m(2);
dT = m(3);
A0 = exp(m(1));

E = [G*m - d];

misfit = E'*diag(wts)*E;











end

