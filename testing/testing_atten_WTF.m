clear all
close all
addpath('/Users/zeilon/Documents/MATLAB/geoff/VpVsQ/JF10Fit')
% parms
Q = 100;
V0 = 5e3; % reference velocity in km/s
rho = 3300;
alpha = 0.27;

Ju = 1./(rho*V0^2);

Temp = 1500;
Z = 10;
gs = 0.001; 

frq = logspace(-2,2,33)';

Vs1 = zeros(size(frq));
Vs2 = zeros(size(frq));
Vs3 = zeros(size(frq));
Qs = zeros(size(frq));

% Calculate V(f) and Qs(f) from JF10 -- experiments
for iff = 1:length(frq)
omega = 2*pi*frq(iff);
[J1,J2] = creep10(Temp,gs,Z/30,omega);

qinv=J2/J1;
gg=1./Ju./sqrt(J1.^2 + J2.^2);

Vs1(iff) = sqrt(gg/rho);
Qs(iff) = 1./qinv;
end

% calculate effective alpha from JF10
G = [log(frq), ones(size(frq))];
d = log(Qs);
m = (G'*G)\G'*d;
alpha_JF10 = m(1)

% compute reference Q and V at 1 second
Q0 = Qs(frq == 1);
V0 = Vs1(frq == 1);
% compute idealised infinite frequency velocity from 1s reference
Vinf = V0*(1 + 1./Q0) 

% Calculate V from physical dispersion relationship for f-independent Q
for iff = 1:length(frq)
Vs2(iff) = V0*( 1 + log(frq(iff))./pi./Q0) ;
end

% Calculate V from physical dispersion relationship for f-dependent Q
% e.g. Lekic 2009 equation 11
for iff = 1:length(frq)
%     qf = (frq^alpha)./Q0;
Vs3(iff) = V0*( 1 + (0.5./Q0)*cot(alpha*pi/2)*(1 - frq(iff).^-alpha) ) ;
end

% from Q0 compute Q if simple f-dependency or constant
Qs2 = Q0*ones(size(frq));
Qs3 = Q0*frq.^alpha;

% compare V predicted by Kanamori 1977 to that coming out of JF10
figure(1), clf, hold on
hv(1) = semilogx(frq,Vs1,'b'); % JF10 velocity
hv(2) = semilogx(frq,Vs2,'r'); % Kanamori f-independent physical dispersion velocity
hv(3) = semilogx(frq,Vs3,'g'); % Anderson f-dependent physical dispersion
set(gca,'xscale','log')
xlabel('freq','Fontsize',12)
ylabel('V','Fontsize',12)
legend(hv,'JF10','Kanamori (a=0)','Anderson (a=0.27)','location','SouthEast')

figure(2), clf, hold on
hq(1) = semilogx(frq,Qs,'b')
hq(2) = semilogx(frq,Qs2,'r')
hq(3) = semilogx(frq,Qs3,'g')
set(gca,'xscale','log','ylim',[20 200])
xlabel('freq','Fontsize',12)
ylabel('Qs','Fontsize',12)
legend(hq,'JF10','Kanamori (a=0)','Anderson (a=0.27)','location','SouthEast')

return
figure(3), clf, hold on
plot(frq,log(frq),'r')
plot(frq,(1-frq.^-alpha)/alpha)
set(gca,'xscale','log')

