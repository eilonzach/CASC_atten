% Script to check my maths and explore some of the consequences of
% frequency dependence for the splitting operator, and inverting it. 
close all
addpath('matguts')
% parms
samprate = 100;

%%  Q & V 
Q0 = 15;
c0 = 4.e3; % reference velocity in km/s

z = 200e3;

alpha = 0.27;

w = 2*pi*logspace(-log10(30),0,20);

y = cot(alpha*pi/2)/2;

Qw = Q0*power(w,alpha);
qinv0 = 1./Q0;

cw = c0*( 1 - (y/Q0)*(w/pi/2).^-alpha);

tstar = z/c0/Q0;

%% Attenuation operator
% Dwt = exp( -0.5*abs(w)*L.*qinv_w./c0 ) .* exp( -1i*w.*L.*(1./c_w - 1./c_inf) );
% Dwt = A exp(-i*w*phi)

A = exp( (-w.*z)./(2.*c0.*Qw) );
A = exp(-tstar.*w.^(1-alpha)/2); %PRECISE

wphi = w.*z.*(1./cw - 1./c0)
% wphi = w.*z.*( (c0*(1 - y.*qinv0.*(w/2/pi).^-alpha)).^-1 - 1/c0)
% wphi = w.*z.*( 1/c0*((1 - y.*qinv0.*(w/2/pi).^-alpha).^-1 - 1))
wphi = w.*z.*( 1/c0*( (1 - y.*qinv0.*(w/2/pi).^-alpha).^-1 - 1 ) )

% % until this point it has been PRECISE, now make Maclurin expansion of
% the (1 - x)^-1  where
x = y.*qinv0.*(w/2/pi).^-alpha;
% wphi_ = w.*z.*( 1/c0*(  1 + x + x.^2 + x.^3 + x.^4          - 1 ) )
% wphi_ = w.*z.*( 1/c0*( x + x.^2 + x.^3 + x.^4 ) )
% wphi_ = w.*z.*( 1/c0*( x + x.^2 + x.^3) )
% wphi_ = w.*z.*( 1/c0*( x + x.^2) )
wphi_ = w.*z.*( 1/c0*( x ) );
% wphi_ = tstar.*y.*(w/2/pi).^-alpha

% wphi_ = tstar.*w.*y.*((2*pi).^alpha).*w.^(-alpha)
% wphi_ = tstar.*w.*y.*((2*pi).^alpha).*w.^(-alpha)
return
%% This approximation rests on how close (1 - x)^-1 is to (1 + x)


%% approximation
Q0s = logspace(1,2.5,10);
for iq = 1:length(Q0s)
    qinv0 = 1./Q0s(iq);
    x = y.*qinv0.*(w/2/pi).^-alpha;
    
    X(iq,:) = x;
    e1(iq,1) = norm( ( (1 - x).^-1) - (1 + x) );
    e2(iq,1) = norm( ( (1 - x).^-1) - (1 + x + x.^2) );
    e3(iq,1) = norm( ( (1 - x).^-1) - (1 + x + x.^2 + x.^3) );
    e4(iq,1) = norm( ( (1 - x).^-1) - (1 + x + x.^2 + x.^3 + x.^4) );

    
    fe1(iq,1) = mean( 1 - (1 + x                     )./( (1 - x).^-1 ) );
    fe2(iq,1) = mean( 1 - (1 + x + x.^2              )./( (1 - x).^-1 ) );
    fe3(iq,1) = mean( 1 - (1 + x + x.^2 + x.^3       )./( (1 - x).^-1 ) );
    fe4(iq,1) = mean( 1 - (1 + x + x.^2 + x.^3 + x.^4)./( (1 - x).^-1 ) );
end

plot(Q0s,[e1,e2,e3,e4])
plot(Q0s,[fe1,fe2,fe3,fe4])

figure(34)
plot(Q0s,X)






