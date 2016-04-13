function [ dtstar_pref,dT_pref,A0_pref,alpha_pref,alpha_misfits,dtstars,dTs,A0s,misfits_amp,misfits_phi ] ...
    = invert_allin1_Aphis_4_STA_dtdtstar_alpha( Amat,phimat,fmids,test_alphas,wtmat,amp2phiwt )
%[ dtstar_pref,dT_pref,A0_pref,alpha_pref,alpha_misfits,dtstars,dTs,A0s, misfits_amp, misfits_phi ] ...
%     = invert_allin1_Aphis_4_STA_dtdtstar_alpha( Amat,phimat,fmids,test_alphas,wtmat,amp2phiwt )
%   Script to simultaneously invert pairwise frequency and phase spectra
%   for dtstar and dT at a whole array of stations, looping over a range of
%   alpha values, solving for the best-fitting alpha and the corresponding
%   dtstar and dT
% 
%  Amat and phimat are matrices each of size [handshake(Nstas) x Nfreq]
%  fmids is a vector of frequencies
%  test_alphas is the vector of alpha values (?0) to test
%  wtmat is a matrix of weights, the same size as Amat
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


%% prelims
if nargin < 5 || isempty(wtmat)
    wtmat = ones(size(Amat));
end
if nargin < 6 || isempty(amp2phiwt)
    amp2phiwt = 1;
end

Npair = size(Amat,1);
Nstas = quadratic_solve(1,-1,-2*Npair);
Nf = length(fmids);
Na = length(test_alphas);

ws = 2*pi*fmids;


%% results structures
dtstars = zeros(Nstas,Na);
dTs = zeros(Nstas,Na);
A0s = zeros(Nstas,Na);

wts = reshape(wtmat',Npair*Nf,1);

d_Amp = reshape(log(Amat)',Npair*Nf,1);
d_phi = reshape(phimat',Npair*Nf,1);

G_Amp = spalloc(Npair*Nf,3*Nstas,4*Nf*Npair);
G_phi = spalloc(Npair*Nf,3*Nstas,4*Nf*Npair);

alpha_misfits = zeros(Na,1);

misfits_amp = zeros(Npair,Na);
misfits_phi = zeros(Npair,Na);

for ia = 1:length(test_alphas)
    alpha = test_alphas(ia);
    
    if alpha==0
        Ax = -0.5*ws;
        Px = 1 - log(ws/2/pi)./pi;
    else
        Ax = -0.5*ws.^(1-alpha);
        Px = 0.5.*(2*pi).^alpha*cot(alpha*pi/2)*ws.^(-alpha);
    end
    
    % data vector: [Npair*Nf x 1] for both A and phi
    % mod parm vector: [A0s(Nstas); dtstar(Nstas); dT(Nstas)]
    
    count = 0;
    for is1 = 1:Nstas
    for is2 = is1+1:Nstas
        count = count+1; % count is the same as the handshake #

        % slot into G matrices
        yind = [1:Nf] + (count-1)*Nf;

        G_Amp(yind,is1)         = -1;
        G_Amp(yind,is2)         = 1;
        G_Amp(yind,Nstas+is1)   = -Ax;
        G_Amp(yind,Nstas+is2)   = Ax;

        G_phi(yind,Nstas+is1)   = -Px;
        G_phi(yind,Nstas+is2)   = Px;
        G_phi(yind,2*Nstas+is1) = -1;
        G_phi(yind,2*Nstas+is2) = 1;
    end 
    end

    constraint_A0  = sparse(1,[1:Nstas]        ,1,1,3*Nstas);
    constraint_Amp = sparse(1,[1:Nstas]+Nstas  ,1,1,3*Nstas);
    constraint_phi = sparse(1,[1:Nstas]+2*Nstas,1,1,3*Nstas);
    
    d_all = [d_Amp;d_phi;0;0;0];
    G_all = [G_Amp;G_phi;constraint_A0;constraint_Amp;constraint_phi];
    
    w_all = [amp2phiwt*wts;wts;1;1;1];
    N = length(w_all);
    spdw_all = sparse(1:N,1:N,w_all,N,N,N);

    
    dw = spdw_all*d_all;
    Gw = spdw_all*G_all;
    
    
%     m = (Gw'*Gw)\Gw'*dw;   
    m = lsqr(Gw,dw,1e-7,100);
    
    
    dtstars(:,ia) = m(Nstas+ [1:Nstas]);
    dTs(:,ia) = m(2*Nstas + [1:Nstas]);
    A0s(:,ia) = exp(m(1:Nstas));

    E = [G_all*m - d_all];
    
    alpha_misfits(ia) = E'*spdw_all*E;
    
    misfits_amp(:,ia) = sum(reshape(wts.*E(           [1:Npair*Nf]).^2,Nf,Npair))';
    misfits_phi(:,ia) = sum(reshape(wts.*E(Npair*Nf + [1:Npair*Nf]).^2,Nf,Npair))';
end

%% minimise misfit
alpha_pref = test_alphas(mindex(alpha_misfits));
dtstar_pref = dtstars(:,mindex(alpha_misfits));
dT_pref = dTs(:,mindex(alpha_misfits));
A0_pref = A0s(:,mindex(alpha_misfits));

figure(77), clf;
plot(test_alphas,alpha_misfits,'-o')
xlabel('F-dependency ($\alpha$)','interpreter','latex','FontSize',22)
ylabel('Global misfit, ($\chi^2$)','interpreter','latex','FontSize',22)
set(gca,'FontSize',14,'box','on')



end

