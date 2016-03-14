function [delta_tstar,delta_T] = combspectra(dat,samprate,parms,ifplot)

if nargin < 4 
    ifplot = false;
end

Tmin = parms.comb.Tmin;
Tmax = parms.comb.Tmax;
Nwds = parms.comb.Nwds;
npol = parms.comb.npol;

pretime = parms.wind.pretime;
prex    = parms.wind.prex;
postx   = parms.wind.postx;
taperx   = parms.wind.taperx;

minacor = parms.qc.minacor;
maxphi  = parms.qc.maxphi;

amp2phiwt = parms.inv.amp2phiwt;
fmax      = parms.inv.fmax;
ifwt      = parms.inv.ifwt;

fnq = samprate/2;
dt = 1./samprate;

Nstas = size(dat,2);
N = handshake(Nstas);

dtstar_all = zeros(N,1);
dT_all = zeros(N,1);

%% prepare filter + cleaning parms
%% Make set of period windows for bandpass filter
Tmids = logspace(log10(Tmin),log10(Tmax),Nwds)';
Twdhs = 0.5*diff(logspace(log10(Tmin/2),log10(2*Tmax),Nwds+1)');
% Twdhs = 0.5*Tmids;
fmids = 1./Tmids;

fprintf('SHOULD PLOT FILTERS IN F-SPACE\n')

%% loop over station pairs
ui = zeros(2*N,1);
uj = zeros(2*N,1);
u  = zeros(2*N,1);
count = 0;

hw = waitbar(count/N,'Progress through station-station dA,dphi');

for is1 = 1:Nstas
for is2 = is1+1:Nstas
    waitbar(count/N,hw)
    count = count+1;

    dat1 = dat(:,is1);
    dat2 = dat(:,is2);

    As = zeros(Nwds,1);
    phis = zeros(Nwds,1);
    wts = zeros(Nwds,1);
    %% loop over frequency bins
    for ii = 1:Nwds
        flo = 1./(Tmids(ii) + Twdhs(ii));
        fhi = 1./(Tmids(ii) - Twdhs(ii));
        fmid = 1./Tmids(ii);
        cp = struct('samprate',samprate,'pretime',pretime,'prex',prex,'postx',postx,...
                    'taperx',taperx,'fhi',fhi,'flo',flo,'npoles',npol,'norm',0);

        [ qdatwf1, qdatf1, qdatwc1, ~, ~, ttws, ~ ] = data_clean( dat1,cp );
        [ qdatwf2, qdatf2, qdatwc2, ~, ~, ~, ~ ] = data_clean( dat2,cp );

        % find observed phase shift
        [dcor, dcstd, dvcstd, acor]=xcortimes([qdatwf1,qdatwf2], dt, pretime, 10,0);
        phi_f_obs = diff(dcor);

        % make phase-corrected time series
        qdatwf2s = interp1(ttws-phi_f_obs,qdatwf2,ttws,'linear',0)';

        % calc. observed amplitude factor
        A_f_obs = (qdatwf2s'*qdatwf1)/(qdatwf1'*qdatwf1);
    %     A_f_obs = 1./((datwf'*qdatwfs)/(datwf'*datwf)); %<< could do it other way too

        % make amplitude-corrected time series
        qdatwf2sa = qdatwf2s./A_f_obs;


        As(ii) = A_f_obs;
        phis(ii) = phi_f_obs;
        if ifwt
            acor = xcorr(qdatwf1,qdatwf2sa,0)^2./(xcorr(qdatwf1,qdatwf1,0)*xcorr(qdatwf2sa,qdatwf2sa,0));
            wts(ii) = acor.^2;
        else
            wts(ii) = 1;
        end
    end

    %% use vectors of amplitude factors and phase shifts to estimate dtstar and dT
    % in theory, A(f)   = exp(pi*f*tstar)
    %  so a plot of ln(A) against freq should have slope pi*dtstar
    %  or a plot of ln(A) against 1/t should have slope pi*dtstar
    % in theory, phi(f) = 1./pi * tstar * ln(fnq/f)
    %  so a plot of phi against ln(freq) should have slope -dtstar./pi and
    %  intercept that related to fnq and dt
    %  or a plot of phi against ln(t) should have slope dtstar./pi

    %% QC
    inds = fmids<=fmax & abs(phis)<maxphi & sqrt(wts)>minacor;

    %% calcs
%     % estimates of dtstar from the data
%     fo1 = fit(fmids(inds),log(As(inds)),'poly1','weight',wts(inds));
%     dtstar_e1 = -fo1.p1./pi;
%     fo2 = fit(log(fmids(inds)),phis(inds),'poly1','weight',wts(inds));
%     dtstar_e2 = -fo2.p1*pi;
%     % estimates of dT from the data
%     dT_anel = fo2.p2 + fo2.p1*log(fnq); % using phase lags, taking into account anelasticity from phase
%     dT_xcor = diff(xcortimes([qdatwc1,qdatwc2], dt, 50, 10,0)); % ignoring anelasticity
% 
%     % estimates from simultaneous inversion of amp and phase data
    [ dtstar,dT,A0,misfit ] = invert_Aphi_4_dtdtstar( As(inds),phis(inds),fmids(inds),fnq, wts(inds),amp2phiwt);


    if ifplot
        figure(3), clf, set(gcf,'pos',[600 10 500,700])
        
        subplot(211), hold on
        scatter(fmids,log(As),70*wts,'or','MarkerFaceColor','r'), set(gca,'Xscale','log')
        plot(fmids,log(A0) - pi*fmids*dtstar,'g'), set(gca,'Xscale','log')
        plot(fmax*[1 1],[-1 1],'--b')
        xlabel('log(freq)','FontSize',18), ylabel('log(Amp)','FontSize',18)
        title(sprintf('Station %.0f vs. station %.0f',is1,is2),'FontSize',20)
        
        subplot(212), hold on
        scatter(fmids,phis,70*wts,'or','MarkerFaceColor','r'), set(gca,'Xscale','log')
        plot(fmids,(log(fnq) - log(fmids))*dtstar./pi + dT,'g'), set(gca,'Xscale','log')
        plot(fmax*[1 1],[-1 1],'--b')        
        xlabel('log(freq)','FontSize',18), ylabel('phase-shift (s)','FontSize',18)
     end
    
    dtstar_all(count) = dtstar;
    dT_all(count) = dT;
    wt_all(count) = sum(inds)./misfit; % weight = 1./misfit, normalised by number of datapoints
%     pause

    ui(2*(count-1)+[1 2]) = count;
    uj(2*(count-1)+[1 2]) = [is1 is2];
    u (2*(count-1)+[1 2]) = [1 -1];
end
end

delete(hw)

G = sparse(ui,uj,u,N,Nstas,2*N);

if ifwt
    W = diag(wt_all); 
else
    W = eye(N);
end

% add constraint
G(N+1,:) = 1;
dtstar_all(N+1,:)=0;
dT_all(N+1,:)=0;
W = diag([diag(W);1]);

delta_tstar = (G'*W*G)\G'*W*dtstar_all;
delta_T = (G'*W*G)\G'*W*dT_all;

cov_dtstar = ( (G'*G)\G' ) * diag(diag(W).^-2) * ( (G'*G)\G' )' ;
std_dtstar = diag(cov_dtstar).^0.5;


end