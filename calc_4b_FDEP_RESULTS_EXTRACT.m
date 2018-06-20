% script to sequentially load the grids of pairwise spectra and solve each
% for dtstar assuming a given alpha (not necessarily zero
% 
% the format of the results file will be the same as for RESULTS_EXTRACT,
%  i.e. a big nstas x nevts matrix of values with nans where unavailable

% clear all
close all
cd /Users/zeilon/Documents/MATLAB/CASC_atten/

ifsave = true;

ALP = 0.5;

phase = 'P';
component = 'Z';

fmids = 1./logspace(log10(1),log10(20),30)';

amp2phiwt = 5;

%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% PAIRS DIRECTORY 
pairdir = '~/Documents/MATLAB/CASC_atten/results_amisfits/'; % needs final slash
% RESULTS DIRECTORY
resdir = '~/Documents/MATLAB/CASC_atten/results/'; % needs final slash


%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

% GET EVENTS+STATIONS DATA
db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,elats,elons,edeps,evtimes,mags] = dbgetv(dbor,'orid','lat','lon','depth','time','ms');
norids = dbnrecs(dbor);
dbsi = dblookup_table(db,'site');
[stas,slats,slons,selevs] = dbgetv(dbsi,'sta','lat','lon','elev');
nstas = dbnrecs(dbsi);
dbclose(db)

all_dtstar_comb = nan(nstas,norids);
all_dT_comb = nan(nstas,norids);
avr = nan(norids,1);
aew = nan(norids,1);

%% re-do the fitting calculations
for ie = 1:norids
    pairfile = sprintf('%s%.0f_amisfits_%s%s.mat',pairdir,orids(ie),phase,component);
    if exist(pairfile,'file')~=2, continue, end
    load(pairfile)
    fprintf('Doing orid %.0f\n',ie)
    
    Amat = pairwise.As;
    phimat = pairwise.phis;
    wtmat = double(pairwise.inds).*pairwise.wts;

%     if ifobs % pull out/use only OBS stations by setting wts to zero for pairs with any non-OBS
%         nstas = length(sts);
%         iob = zeros(nstas,1);
%         hob = zeros(handshake(nstas),1);
%         for is = 1:nstas, iob(is) = ~isempty(which_OBS(sts{is})); end
%         k=0;
%         for is1 = 1:nstas
%         for is2 = is1+1:nstas
%         k=k+1;
%         wtmat(k,:) = wtmat(k,:)*double(iob(is1) & iob(is2));
%         end
%         end
%     end

    %% Calculate dtstar etc.
    [ delta_tstar_pref,delta_T_pref,A0_pref,alpha_pref,alpha_VR,alpha_Ew ] ...
        = calc_fdependent( Amat,phimat,fmids,ALP,wtmat,amp2phiwt,2,'',0);

    %% sort into results matrix
    [~,ind,~] = intersect(stas,sts,'stable');
    all_dtstar_comb(ind,ie) = delta_tstar_pref;
    all_dT_comb(ind,ie) = delta_T_pref;
    avr(ie) = alpha_VR;
    aew(ie) = alpha_Ew;
    
end

fprintf('Mean variance reduction is %.3f\n',nanmean(avr));
fprintf('Summed weighted error is %.2f\n',nansum(aew));


if ifsave
    resfile_dtstar = sprintf('all_dtstarcomb_%s_%s_alp%03.0f',phase,component,100*ALP);
    resfile_dT = sprintf('all_dTcomb_%s_%s_alp%03.0f',phase,component,100*ALP);

    save([resdir,resfile_dtstar],'all_dtstar_comb')
    save([resdir,resfile_dT],'all_dT_comb')
end
return
