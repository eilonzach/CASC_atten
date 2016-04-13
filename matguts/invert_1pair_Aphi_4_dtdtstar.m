function [ dtstar,dT,A0,misfit, E, misfit_amps,misfit_phis ] = invert_1pair_Aphi_4_dtdtstar( As,phis,freqs,wts,amp2phiwt,alpha)
% [ dtstar,dT,A0,misfit, E, misfit_amps,misfit_phis ] = invert_1pair_Aphi_4_dtdtstar( As,phis,freqs,[wts=1],[amp2phiwt=1],[alpha=0])
%   Script to invert frequency and phase spectra for dtstar and dT
% 
% if Q is frequency independent (alpha==0)
%     A = A0*exp(-pi*dtstar*f)
%     ln(A) = ln(A0) - pi*dtstar*f 
%     phi = (ln(f) - ln(fNq))*dtstar/pi + dT
% 
% elseif Q is frequency dependent (alpha~   =0)
%     A = A0*exp(-(pi/((2*pi)^alpha)) * f^(1-alpha) * dtstar)
%     ln(A) = ln(A0) - (pi/((2*pi)^alpha)) * f^(1-alpha) * dtstar
%     phi = 0.5*cot(alpha*pi/2)*f^alpha + dT

if nargin < 4 || isempty(wts)
    wts = ones(size(As));
end
if nargin < 5 || isempty(amp2phiwt)
    amp2phiwt = 1;
end
if nargin < 6
    alpha = 0;
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
if alpha==0
    G([1:N],2) = -0.5*ws;
    G(N+[1:N],2) = 1 - log(ws/2/pi)./pi;
elseif alpha~=0
%     G([1:N],2) = -0.5*( ws.^(1-alpha) );
%     G(N+[1:N],2) = ws.^(-alpha) .* (2*pi).^(1-alpha) * cot(alpha*pi/2);
    G([1:N],2) = -0.5*ws.^(1-alpha) ;
    G(N+[1:N],2) = 0.5*ws.^(-alpha) .* (2*pi).^alpha * cot(alpha*pi/2);
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

