% script to sequentially load the grids of pairwise spectra and solve each
% for dtstar assuming a given alpha (not necessarily zero
% 
% the format of the results file will be the same as for RESULTS_EXTRACT,
%  i.e. a big nstas x nevts matrix of values with nans where unavailable

% clear all
% close all
cd /Users/zeilon/Documents/MATLAB/CASC_atten/

ifsave = true;

alp = 0.27;

alps = [0.27]';
Ealp = zeros(size(alps));
for ia = 1:length(alps);
alp = alps(ia);

phacomp = 'PZ';

fmids = 1./logspace(log10(1),log10(20),30)';
Nf = length(fmids);

amp2phiwt = 5;

%% QC parms
fmax = 1/3;
fwt = double(fmids<=fmax);
maxdist = 5; % in degrees
minWtsum = 5;
ifonlyOBS = 0.5; % or do 0.5 to allow if one of two is obs
oridfirst = 50;
oridlast = 300;

parms = struct('alp',alp,'amp2phiwt',amp2phiwt,'fmax',fmax,'maxdist_deg',maxdist,...
               'minWtsum',minWtsum,'ifonlyOBS',ifonlyOBS,'oridfirst',oridfirst,'oridlast',oridlast);

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

%% starts
ws = 2*pi*fmids;
w0 = 2*pi;
if alp==0
    Ax = -0.5*ws;
    Px = 1 - log(ws/2/pi)./pi;
else
    Ax = -0.5 * w0^alp * ws.^(1-alp);
    Px = 0.5 * cot(alp*pi/2) * (ws/w0).^(-alp);
end

% GET EVENTS+STATIONS DATA
[ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam );
[ Nstas,stas,slats,slons,selevs,ondate,offdate,staname,statype,network ] = db_stadata( dbdir,dbnam );
sages = jdf_crust_age( slats,slons );

%% distances
% for is1 = 1:Nstas
% for is2 = is1+1:Nstas
%     dd(is1,is2) = distance(slats(is1),slons(is1),slats(is2),slons(is2));
% end
% end
% dd = [dd;zeros(1,Nstas)];
% dd = dd+dd';
% save('mapdata/sta_sta_dist.mat','dd');
load('mapdata/sta_sta_dist.mat'); % distance in degrees


%% results structures
% max N of measurements (this is the number of pairs, not of amplitudes/phases measured
Npair_max = handshake(round(Nstas/10))*norids; % (say on average 1/10 of stas are operational at any time - if needed, matrices will lengthen.

d_Amp = nan(Nf*Npair_max,1);
d_phi = nan(Nf*Npair_max,1);
datinfo = nan(Nf*Npair_max,3);
Wts  = nan(Nf*Npair_max,1); % inverse of the weights

% N will be Npair(=count)*Nf 
% M will be 3*Nstas for lnAmp, dtstar, dT (in that order)
% so size of G will be [Npair*Nf x 3*Nstas] 

G_Amp = spalloc(Nf*Npair_max,3*Nstas,4*Nf*Npair_max);
G_phi = spalloc(Nf*Npair_max,3*Nstas,4*Nf*Npair_max);

count = 0; % the count is the number of station pairs ==> length of G will be Nf*count
for ie = oridfirst:oridlast
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
        
        % QCs
        if dd(sind1,sind2) > maxdist % don't include pair if too far apart
            continue
        end
        if sum(wtmat(evcount,:)' .* fwt) < minWtsum, fprintf('badwt\n'),continue, end
        if ifonlyOBS==1
            if any(~strcmp(statype([sind1,sind2]),'OBS')), continue, end
        elseif ifonlyOBS==0.5
            if all(~strcmp(statype([sind1,sind2]),'OBS')), continue, end
        end

        count = count+1;       % overall count
        
        yind = [1:Nf] + (count-1)*Nf;
        
        % add elements of G matrices
        G_Amp(yind,sind1)         = -1;
        G_Amp(yind,sind2)         = 1;
        G_Amp(yind,Nstas+sind1)   = -Ax;
        G_Amp(yind,Nstas+sind2)   = Ax;

        G_phi(yind,Nstas+sind1)   = -Px;
        G_phi(yind,Nstas+sind2)   = Px;
        G_phi(yind,2*Nstas+sind1) = -1;
        G_phi(yind,2*Nstas+sind2) = 1;

        d_Amp(yind) = log(Amat(evcount,:))';
        d_phi(yind) = phimat(evcount,:)';
        Wts(yind)  = wtmat(evcount,:)' .* fwt;
        
        datinfo(yind,1) = sind1; datinfo(yind,2) = sind2; datinfo(yind,3) = ie;
    end % loop on sta2
    end % loop on sta1
end % loop on orids
%% delete nans
G_Amp = G_Amp(1:Nf*count,:);
G_phi = G_phi(1:Nf*count,:);
d_Amp = d_Amp(1:Nf*count,:);
d_phi = d_phi(1:Nf*count,:);
Wts  = Wts(1:Nf*count,:);
datinfo = datinfo(1:Nf*count,:);

%% find stations with data
fprintf('Stripping out nans\n')
isdat_Amp = reshape(full(any(G_Amp))',Nstas,3)';
isdat_phi = reshape(full(any(G_phi))',Nstas,3)';
gdstas = any([isdat_Amp;isdat_phi])';


%% correct for repeat stations
stinds = find(gdstas);
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
            if dd(sind1,sind2) < 0.1 % distance is less than 0.1 deg
            % we have a match!
            doubles = [doubles;sind1,sind2];
            done = [done;sind2];
            end
        end
        
    end
end
% correct doubles stas in Gs
for is = 1:length(doubles)
    xind1 = doubles(is,1) + [0,Nstas,2*Nstas];
    xind2 = doubles(is,2) + [0,Nstas,2*Nstas];
    G_Amp(:,xind1) = G_Amp(:,xind1) + G_Amp(:,xind2);
    G_phi(:,xind1) = G_phi(:,xind1) + G_phi(:,xind2);
end
% wipe doubled stas in gdstas
gdstas(done) = false;
stinds = find(gdstas);

%% subset Gs to stations with gdstas (i.e. there are obs + not doubled)
Ngd = sum(gdstas);
xind = find(gdstas)*[1 1 1] + ones(Ngd,1)*[0,Nstas,2*Nstas];
xind = xind(:)'; % these are the indices of all columns in Gs to keep

G_Amp_do = G_Amp(:,xind);
G_phi_do = G_phi(:,xind);

Nobs = full(sum(G_Amp_do(:,1:Ngd)~=0))'/Nf; % count N of obs per station, summing non-zeros down rows and dividing by Nf fmids
Nwobs = full(sum(sparse(1:Nf*count,1:Nf*count,Wts,Nf*count,Nf*count,Nf*count)*G_Amp_do(:,1:Ngd)~=0))'/Nf; % count weighted obs per station, multiplying by weights (some of which are zero) and then summing non-zeros down rows and dividing by Nf fmids

%% make final matrices
fprintf('Making matrices\n')
constraint_lnAmp  = sparse(1,[1:Ngd]      ,1,1,3*Ngd);
constraint_dtstar = sparse(1,[1:Ngd]+Ngd  ,1,1,3*Ngd);
constraint_dT     = sparse(1,[1:Ngd]+2*Ngd,1,1,3*Ngd);

d_all = [d_Amp;d_phi;0;0;0];
G_all = [G_Amp_do;G_phi_do;constraint_dtstar;constraint_dT;constraint_lnAmp];

w_all = [amp2phiwt*Wts;Wts;1;1;1];
N = length(w_all);
spdw_all = sparse(1:N,1:N,w_all,N,N,N);


dw = spdw_all*d_all;
Gw = spdw_all*G_all;

m = lsqr(Gw,dw,1e-6,1000);
lnA = m(1:Ngd);
dtstar = m([1:Ngd] + Ngd);
dT = m([1:Ngd] + 2*Ngd); 


%% kill useless rows for the bootstrap
killrow = find((sum(Gw.^2,2)==0) & (dw==0));
Gw(killrow,:) = [];
dw(killrow,:) = [];
Ndat = length(dw);

%% Start bootstrap!
Nboot = 1000;
lnAs = zeros(Ngd,Nboot);
dtstars = zeros(Ngd,Nboot);
dTs = zeros(Ngd,Nboot);
return
% % indices from balanced resampling
% inds = repmat([1:Ndat]',Nboot,1);
% error('DO NOT DO THIS! TAKES TOO LONG AND KILLS MATLAB')
% rp = randperm(length(inds));
% inds = inds(rp);
% inds = reshape(inds,Nboot,Ndat);
% 

fprintf('Bootstrapping errors for slope fit\n')
for ib = 1:Nboot
    ib
indx = unidrnd(Ndat,Ndat,1);
% solve lsqr problem
m = lsqr(Gw(indx,:),dw(indx),1e-6,1000);

% parse results
lnAs(:,ib) = m(1:Ngd);
dtstars(:,ib) = m([1:Ngd] + Ngd);
dTs(:,ib) = m([1:Ngd] + 2*Ngd); 

end

% Find 2-sigma upper and lower bounds 
dtstars_ub = zeros(Ngd,1);
dtstars_lb = zeros(Ngd,1);
for is = 1:Ngd
    dts_sort = sort(dtstars(is,:));
    dtstars_ub(is) = dts_sort(round(Nboot*0.95));
    dtstars_lb(is) = dts_sort(round(Nboot*0.05));
end
    


% mopts_sort = sort(mopts);bopts_sort = sort(bopts);
% mup = mopts_sort(round(Nboot*0.95)) %4.96
% mlo = mopts_sort(round(Nboot*0.05)) %2.97
% bup = bopts_sort(round(Nboot*0.95)) %0.00
% blo = bopts_sort(round(Nboot*0.05)) %-0.30

%% analyse error

ew = (Gw*m - dw);
Ealp(ia) = ew'*ew;

end % loop on alphas

plot(alps,Ealp,'o-')
return
%% simple plot
figure(1); clf, set(gcf,'pos',[200 200 600 800])
latlims=[39 51]; lonlims=[-132 -120];
subplot(4,1,1:3), hold on
scatter(slons(stinds),slats(stinds),Nwobs+0.001,dtstar,'filled')
colormap(parula), caxis([-2 2])
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
subplot(4,1,4), hold on
ind = Nwobs>40 & selevs(stinds)<-1000;
scatter(sages(stinds(ind)),dtstar(ind),Nwobs(ind)/2+0.001,'r','filled')
grid on
if any(doubles),
    [~,ind,~] = intersect(stinds,doubles(:,1),'stable');
    scatter(sages(stinds(ind)),dtstar(ind),Nwobs(ind)/2+0.001,'b')  
end
xlim([0 11]), ylim([-3 2])

%% collate results
results = struct('dtstar',dtstar,'dtstar_lb',dtstars_lb,'dtstar_ub',dtstars_ub,...
                 'dT',dT,'lnA',lnA,'Nobs',Nobs,'Nwobs',Nwobs,...
                 'stas',{stas(stinds)},'slats',slats(stinds),'slons',slons(stinds),...
                 'selevs',selevs(stinds),'sages',sages(stinds),'statypes',{statype(stinds)},...
                 'parms',parms);

if ifsave
    resfile = sprintf('allin1_stav_dtstar_wbounds_%s_alp%03.0f',phacomp,100*alp);
    fprintf('Saving!\n')
    save([resdir,resfile],'results')
end

return


%% check some fits

% plot...
gdevts = unique(datinfo(:,3));
figure(44), clf, set(gcf,'pos',[ 440 -139 1430 797])
figure(45), clf, set(gcf,'pos',[ 440 -139 1430 797])
kk = 0;
while kk < 7*8
    % pick random pair of stations   
    rs = randsample(Ngd,2); stind1 = stinds(rs(1)); stind2 = stinds(rs(2));
    % pick random event
    re = randsample(length(gdevts),1); ie = gdevts(re);
    % find rows of data matrices where these are
    ii = find( min(datinfo(:,[1,2]),[],2)==min(stind1,stind2) & max(datinfo(:,[1,2]),[],2)==max(stind1,stind2) & datinfo(:,3)==ie);
    if isempty(ii), continue, end
    % grab differential values according to our results

    [Amat_pred,phimat_pred] = pred_Amat_phimat( dtstar(rs),dT(rs),10.^(lnA(rs)),fmids,alp );
    
    kk = kk + 1;
    
    figure(44), subplot(8,7,kk), hold on, 
    scatter(fmids,d_Amp(ii),50*Wts(ii) + 0.01,'or','filled'),
    plot(fmids,log(Amat_pred),'Linewidth',2), 
    title(sprintf('Stas %.0f - %.0f, ev%.0f',rs,re))

    figure(45), subplot(8,7,kk), hold on, 
    scatter(fmids,d_phi(ii),50*Wts(ii) + 0.01,'og','filled'),
    plot(fmids,phimat_pred,'Linewidth',2), 
    title(sprintf('Stas %.0f - %.0f, ev%.0f',rs,re))

end



