% script to sequentially load the grids of pairwise spectra and solve each
% for dtstar assuming a given alpha (not necessarily zero
% 
% the format of the results file will be the same as for RESULTS_EXTRACT,
%  i.e. a big nstas x nevts matrix of values with nans where unavailable

clear all
close all
cd /Users/zeilon/Documents/MATLAB/CASC_atten/

ifsave = false;

alp = 0.27;
phacomp = 'ST';

fmids = 1./logspace(log10(1),log10(20),30)';
fmax = 1/3;

amp2phiwt = 5;
dampwt = 0.00001;
maxdist = 100; % in degrees


%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% PAIRS DIRECTORY 
pairdir = '~/Documents/MATLAB/CASC_atten/results_pairspecs/'; % needs final slash
% RESULTS DIRECTORY
resdir = '~/Documents/MATLAB/CASC_atten/results/'; % needs final slash


%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

% GET EVENTS+STATIONS DATA
[ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam );
[ nstas,stas,slats,slons,selevs,ondate,offdate,staname,statype,network ] = db_stadata( dbdir,dbnam );
sages = jdf_crust_age( slats,slons );

% max N of measurements
Nmax = handshake(round(nstas/4))*norids; % (say at most 1/4 of stas are operational at any time

% prep G mat elements
ui = nan(2*Nmax,1);
uj = nan(2*Nmax,1);
u  = nan(2*Nmax,1);
count = 0;

dtstar_pairwise = nan(Nmax,1);
dT_pairwise = nan(Nmax,1);
lgA0_pairwise = nan(Nmax,1);
Winv = nan(Nmax,1); % inverse of the weights

for ie = 50:300
    pairfile = sprintf('%s%.0f_pairspecs_%s.mat',pairdir,orids(ie),phacomp);
    if exist(pairfile,'file')~=2, continue, end
    load(pairfile)
    fprintf('Doing orid %.0f\n',ie)
	
    Amat = pairwise.As;
    phimat = pairwise.phis;
    wtmat = double(pairwise.inds).*pairwise.wts;
    
    evNstas = unhandshake(size(Amat,1));
    [~,evINDsta] = intersect(stas,sts,'stable');

    evcount = 0;
    for is1 = 1:evNstas
    for is2 = is1+1:evNstas
        evcount = evcount+1;   % handshake # within this event
        
        sind1 = evINDsta(is1); % overall index of station 1
        sind2 = evINDsta(is2); % overall index of station 2
        
        if distance(slats(sind1),slons(sind1),slats(sind2),slons(sind2)) > maxdist % don't include pair if too far apart
            continue
        end
        
        [ dtstar,dT,A0,misfit] ...
            = invert_1pair_Aphi_4_dtdtstar(Amat(evcount,:),phimat(evcount,:),fmids,...
                        double(fmids<=fmax)'.*wtmat(evcount,:),amp2phiwt,alp);

        if isnan(dtstar)
            continue
        end
        
        count = count+1;       % overall count
        
        % make elements of eventual G matrix
        ui(2*(count-1)+[1 2]) = count;
        uj(2*(count-1)+[1 2]) = [sind1 sind2];
        u (2*(count-1)+[1 2]) = [-1 1]; % delta is value of 2 - value of 1 + event term
        
        dtstar_pairwise(count) = dtstar;
        dT_pairwise(count) = dT;
        lgA0_pairwise(count) = log(A0);
        Winv(count) = misfit./sum(wtmat(evcount,:)~=0); % weight  will be  1./misfit, normalised by number of datapoints
    end 
    end
end

%% delete nans
fprintf('Stripping out nans\n')
u(2*count+1:end) = [];
ui(2*count+1:end) = [];
uj(2*count+1:end) = [];

dtstar_pairwise(count+1:end) = [];
dT_pairwise(count+1:end) = [];
lgA0_pairwise(count+1:end) = [];
Winv(count+1:end) = [];

d1 = dtstar_pairwise;
d2 = dT_pairwise;
d3 = lgA0_pairwise;


% correct for repeat stations/ none stations
stinds = unique(uj);
Ngd = length(stinds);
doubles = [];
done = [];
for is1 = 1:Ngd
    sind1 = stinds(is1);
    if any(done==sind1), continue; end
    for is2 = is1+1:Ngd
        sind2 = stinds(is2);
        if any(done==sind2), continue; end
        namel = length(stas{sind1}); if namel<4, continue; end % don't do for short sta names
        if strcmp(stas{sind1}(1:end-1),stas{sind2}(1:end-1)) %  1:N-1 characters in name match
            if distance(slats(sind1),slons(sind1),slats(sind2),slons(sind2)) < 0.1 % distance is less than 0.1 deg
            % we have a match!
%             fprintf('%s %s\n',stas_parsed{is1},stas_parsed{is2})
            doubles = [doubles;sind1,sind2];
            done = [done;sind2];
            end
        end
    end
end
% replace doubles stas with this
for is = 1:length(doubles)
    uj(uj==doubles(is,2)) = doubles(is,1);
end
% count Nobs/sta
Nobs = zeros(length(stas),1);
for is = 1:nstas
    Nobs(is) = sum(uj == is);
end 
    
%% solve the least squares problem
fprintf('Making matrices\n')

G = sparse(ui,uj,u,count,nstas,2*count);

W = sqrt(1./Winv); 

% damping - need for underdetermined stations
L = sparse([1:nstas],[1:nstas],ones(nstas,1),nstas,nstas,nstas);

% constraint line - station terms must average to zero;
C = sparse([ones(1,nstas)]);

F = [G; dampwt*L; 1e-9*C];
f1 = [d1; zeros(nstas,1); 0];
f2 = [d2; zeros(nstas,1); 0];
f3 = [d3; zeros(nstas,1); 0];

Nmax = size(F,1);
Wmat = sparse([1:Nmax],[1:Nmax],[W; ones(nstas,1); 1],Nmax,Nmax,Nmax);

W2 = 1./Winv; 
Wmat = sparse([1:Nmax],[1:Nmax],[W2; ones(nstas,1); 1],Nmax,Nmax,Nmax);
% m1 = [F'*W2*F]\F'*W2*f1;


fprintf('Solving weighted least squares\n')

Fw = Wmat*F; 
f1w = Wmat*f1;

m1 = lsqr(F,f1,1e-6,1000);
% m1 = lsqr(Fw,f1w,1e-6,1000)

    
m1(abs(m1)<0.01) = nan;

figure(1); clf, set(gcf,'pos',[200 200 600 800])
latlims=[39 51]; lonlims=[-132 -120];
subplot(4,1,1:3), hold on
scatter(slons,slats,Nobs/2+0.001,m1,'filled')
colormap(parula), caxis([-2 1])
% [chgrdX,chgrdY] = meshgrid(linspace(lonlims(1),-123,50),linspace(40,latlims(2),60)); 
% [ chgrdage,chrons,F ] = jdf_crust_age(chgrdY,chgrdX);
% % chtick = unique(chrons.age(~isnan(chrons.age))); % plot chrons
% chtick = 1:12; % plot Ma
% for ic = 1:length(chtick)
%     ch = chtick(ic);
%     contour(chgrdX,chgrdY,chgrdage,[ch ch],'LineWidth',2,'linestyle','--',...
%         'color',colour_get(ch,max(chtick),min(chtick),hot))
% end
xlim([-132 -120]), ylim(latlims)


% caxis([-0.8 0.5])
subplot(4,1,4)
scatter(sages(Nobs>100),m1(Nobs>100),Nobs(Nobs>100)/2+0.001,'r','filled')



return

if ifsave
    resfile_dtstar = sprintf('all_dtstarcomb_%s_%s_alp%03.0f',phase,component,100*alp);
    resfile_dT = sprintf('all_dTcomb_%s_%s_alp%03.0f',phase,component,100*alp);

    save([resdir,resfile_dtstar],'all_dtstar_comb')
    save([resdir,resfile_dT],'all_dT_comb')
end
return
