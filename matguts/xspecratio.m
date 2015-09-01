% function [ delta_tstar ] = xspecratio( specs,frq,fmax,wtopt )
%[ delta_tstar ] = Untitled( specs,frq,fmax,wtopt )
%   Function to calculate differential t-star by taking pair-wise spectral
%   ratios for all combinations of stations' spectra and then using a
%   least squares approach to solve for best-fitting delta t-star with the
%   constraint that the mean dtstar=0. 
%   N.B. this script assumes that the spectral ratio can be fit by a linear
%   relationship, and the fit is over the frequency interval 0 to fmax
% 
%  INPUTS:
%    specs   - [nfreq x nstas] matrix of all the spectra, in columns
%    frq     - [nfreq x 1] vector of frequencies 
%    fmax    - high frequency end of fitting window
%    wtopt   - (optional) option to weight least-squares inversion by
%               formal error in the linear fit to the spectral ratio
% 
%  OUTPUT:
%    delta_tstar - 

hifrq = .1;
plotopt = 1;
figure(2), clf, hold on

% if nargin < 4
%     wtopt = 0;
% end

M = size(specss,2);
N = handshake(M);

dtst = zeros(N,1);
wtst = zeros(N,1);

ui = zeros(2*N,1);
uj = zeros(2*N,1);
u  = zeros(2*N,1);
count = 0;

for ii = 1:M
for jj = ii+1:M
    count = count+1;
    % spectral ratio
    R = specss(:,ii)./specss(:,jj);
    % take logarithm
    lnR = log(R); 
    ind = frq < hifrq;


    
    fo = fit(frq(ind),lnR(ind),'poly1');
    cint = confint(fo)/pi; 
    
    dtst(count) = fo.p1/(-pi);
    wtst(count) = diff(cint(1,:)'); 
    
    if plotopt
    hold on
    hp = plot(frq(ind),lnR(ind),'o-');
	set(hp,'color',colour_get(eqar(indgd(ii)).slon,max([eqar.slon]),min([eqar.slon])))
    xlim([0 hifrq])
%     plot(fo)
    end
    
    ui(2*(count-1)+[1 2]) = count;
    uj(2*(count-1)+[1 2]) = [ii jj];
    u (2*(count-1)+[1 2]) = [1 -1];
end
end
G = sparse(ui,uj,u,N,M,2*N);

if wtopt
    W = diag(wtst.^-1);
else
    W = eye(N);
end

% add constraint
G(N+1,:) = 1;
dtst(N+1)=0;
W = diag([diag(W);1]);

delta_tstar = (G'*W*G)\G'*W'*dtst;

% end

