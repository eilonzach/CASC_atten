function [ dtstar,dT,A0,misfit, E, misfit_amps,misfit_phis ] = invert_1pair_Aphi_4_dtdtstar( As,phis,freqs,wts,amp2phiwt,alp,w0)
% [ dtstar,dT,A0,misfit, E, misfit_amps,misfit_phis ] = invert_1pair_Aphi_4_dtdtstar( As,phis,freqs,[wts=1],[amp2phiwt=1],[alp=0],[w0=2*pi])
%   Script to invert frequency and phase spectra for dtstar and dT
% 
% if Q is frequency independent (alp==0)
%     A = A0*exp(-pi*dtstar*f)
%     ln(A) = ln(A0) - pi*dtstar*f 
%     phi = (ln(f) - ln(fNq))*dtstar/pi + dT
% 
% elseif Q is frequency dependent (alp~   =0)
%     A = A0*exp(-pi*f0^alp * f^(1-alp) * dtstar)
%     ln(A) = ln(A0) - pi*f0^alp * f^(1-alp) * dtstar
%     phi = 0.5*cot(alp*pi/2)*(f/f0)^alp + dT

if nargin < 4 || isempty(wts)
    wts = ones(size(As));
end
if nargin < 5 || isempty(amp2phiwt)
    amp2phiwt = 1;
end
if nargin < 6
    alp = 0;
end
if nargin < 7
    w0 = 2*pi;
end

As = As(:); 
phis = phis(:); 
wts = wts(:);

ws = 2*pi*freqs;

N = length(As);

% Model parameter vector will be [ln(A0); dtstar; dT]

% data vector - amplitudes and phis
d = [log(As);phis];
% G matrix - relationship between amplitudes, phis, and mod parms
G = zeros(2*N,3);
G([1:N],1) = 1;
if alp==0
    G([1:N],2) = -0.5*ws;
    G(N+[1:N],2) = 1 - log(ws/2/pi)./pi;
elseif alp~=0
%     G([1:N],2) = -0.5*( ws.^(1-alp) );
%     G(N+[1:N],2) = ws.^(-alp) .* (2*pi).^(1-alp) * cot(alp*pi/2);
    G([1:N],2) = -0.5*ws.^(1-alp) ;
    G(N+[1:N],2) = 0.5*(ws/w0).^(-alp) * cot(alp*pi/2);
end
G(N+[1:N],3) = 1;
    
% weight
dw = diag([amp2phiwt*wts;wts])*d;
Gw = diag([amp2phiwt*wts;wts])*G;

m = (Gw'*Gw)\Gw'*dw;    

dtstar = m(2);
dT = m(3);
A0 = exp(m(1));

lnApred = G(1:N,:)*m;
Apred = exp(lnApred);

phipred = G(1+N:2*N,:)*m;

Ea = Apred - As;        % amplitudes error (in abs amp space, not log!
Ep = phipred - phis;    % phis error

E = [Ea;Ep]; % error vector

misfit = E'*diag([wts;wts])*E; % report weighted misfit (no amp2phiwt)

misfit_amps = Ea'*diag(wts)*Ea;
misfit_phis = Ep'*diag(wts)*Ep;











end

