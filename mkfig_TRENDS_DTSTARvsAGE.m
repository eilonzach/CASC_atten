% cycle through events and calculate differential travel times for all stations
% clear all
close all
%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% RESULTS DIRECTORY
resdir = '~/Documents/MATLAB/CASC_atten/results/'; % needs final slash
cd ~/Documents/MATLAB/CASC_atten/
addpath('matguts')

%% parameters
ifsave = 1;
ifpreds = 0;
method = 'comb';% 'comb' or 'specR'
tstlim = 1.*[-1.0 0.5];

km_per_ma = 40; % oldest subducting JdF crust is 12 Ma old, 480 km from the ridge


%% compute synthetic tstar
if ifpreds
    run('synthetic_tstar/pred_tstar'), close all
end
% run('synthetic_tstar/pred_tstar_WET')

%% GET EVENTS+STATIONS DATA
[ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam );
[ nstas,stas,slats,slons,selevs,~,~,~,statype ] = db_stadata( dbdir,dbnam );

% Work out distance to ridge, and age
[~,Xrdg,gdjdf] = dist_to_rdg(slats,slons);
age = Xrdg/km_per_ma;


%% LOAD RESULTS
load([resdir,'stav_dTxcorr_P_Z.mat']); dTP = results; 
load([resdir,'stav_dTxcorr_S_T.mat']); dTS = results;
load([resdir,'stav_dtstar',method,'_P_Z.mat']); dtsP = results;
load([resdir,'stav_dtstar',method,'_S_T.mat']); dtsS = results;
% find which obs
for is = 1:length(dTP.stas), dTP.isob(is,1) = strcmp('OBS',statype(strcmp(stas,strtok(dTP.stas(is),'_')))); end
for is = 1:length(dTS.stas), dTS.isob(is,1) = strcmp('OBS',statype(strcmp(stas,strtok(dTS.stas(is),'_')))); end
for is = 1:length(dtsP.stas), dtsP.isob(is,1) = strcmp('OBS',statype(strcmp(stas,strtok(dtsP.stas(is),'_')))); end
for is = 1:length(dtsS.stas), dtsS.isob(is,1) = strcmp('OBS',statype(strcmp(stas,strtok(dtsS.stas(is),'_')))); end
% calc. age from distance to ridge
[~,dTP.Xrdg,dTP.ijdf] = dist_to_rdg(dTP.slats,dTP.slons);
[~,dTS.Xrdg,dTS.ijdf] = dist_to_rdg(dTS.slats,dTS.slons);
[~,dtsP.Xrdg,dtsP.ijdf] = dist_to_rdg(dtsP.slats,dtsP.slons);
[~,dtsS.Xrdg,dtsS.ijdf] = dist_to_rdg(dtsS.slats,dtsS.slons);


%% ==================================================================== %% 
%% =============================  DTSTAR  ============================= %% 
%% ==================================================================== %% 


% % LOAD RESULTS
% if strcmp(method,'specR')
%     load([resdir,'all_dtstar_P_Z.mat']); Pall_dtstar = all_dtstar;
%     load([resdir,'all_dtstar_S_T.mat']); Sall_dtstar = all_dtstar;
% elseif strcmp(method,'comb')
%     load([resdir,'all_dtstarcomb_P_Z.mat']); Pall_dtstar = all_dtstar_comb;
%     load([resdir,'all_dtstarcomb_S_T.mat']); Sall_dtstar = all_dtstar_comb;
% end
% 
% [ Psta_dtstar] = lsq_sta_evt( Pall_dtstar,0.01);
% Pstd = nanvariance(Pall_dtstar,false,2).^(1/2);
% PNobs = sum(~isnan(Pall_dtstar),2);
% 
% [ Ssta_dtstar] = lsq_sta_evt( Sall_dtstar,0.01);
% Sstd = nanvariance(Sall_dtstar,false,2).^(1/2);
% SNobs = sum(~isnan(Sall_dtstar),2);

%% Work out good stations, and plot them
%P
Pgd  = dtsP.dtstar~=0 & abs(dtsP.dtstar)<3 & dtsP.Nobs>=5;
gdjdf  = Pgd & dtsP.isob & dtsP.ijdf;
gdgor  = Pgd & dtsP.isob & ~dtsP.ijdf;
dcshift = 0.4;
figure(31),clf, set(gcf,'pos',[300 500 800 360]),hold on
if ifpreds
for ig = 1:length(gs)
hp(ig) = plot(ages,dtstar_P(:,ig)-min(dtstar_P(:,ig)),'k','LineWidth',2,'LineStyle',gsls{ig});
end
end
scatter(abs(dtsP.Xrdg(gdgor)/km_per_ma),dtsP.dtstar(gdgor)+dcshift,20*dtsP.Nobs(gdgor),...
    'k','MarkerFaceColor','b','LineWidth',1.5)
scatter(abs(dtsP.Xrdg(gdjdf)/km_per_ma),dtsP.dtstar(gdjdf)+dcshift,20*dtsP.Nobs(gdjdf),...
    'k','MarkerFaceColor','r','LineWidth',1.5)
if ifpreds
    legend(hp,'0.1mm','1mm','1cm','location','northeast')
end

xlabel('\textbf{Age of oceanic lithosphere (Ma)}','Interpreter','latex','FontSize',22)
ylabel('$\mathbf{\Delta t^{\ast}_P}$','Interpreter','latex','FontSize',22)
set(gca,'ylim',[-0.5 1.0],'xlim',[0 12],'fontsize',16,'box','on','linewidth',2)

%S
Sgd  = dtsS.dtstar~=0 & abs(dtsS.dtstar)<3 & dtsS.Nobs>=5;
gdjdf  = Sgd & dtsS.isob & dtsS.ijdf;
gdgor  = Sgd & dtsS.isob & ~dtsS.ijdf;
dcshift = 2.0;
figure(32),clf, set(gcf,'pos',[300 000 800 360]),hold on
if ifpreds
for ig = 1:length(gs)
hs(ig) = plot(ages,dtstar_S(:,ig)-min(dtstar_S(:,ig)),'k','LineWidth',2,'LineStyle',gsls{ig});
end
end
scatter(abs(dtsS.Xrdg(gdgor)/km_per_ma),dtsS.dtstar(gdgor)+dcshift,20*dtsS.Nobs(gdgor),...
    'k','MarkerFaceColor','b','LineWidth',1.5)
scatter(abs(dtsS.Xrdg(gdjdf)/km_per_ma),dtsS.dtstar(gdjdf)+dcshift,20*dtsS.Nobs(gdjdf),...
    'k','MarkerFaceColor','r','LineWidth',1.5)
if ifpreds
    legend(hs,'0.1mm','1mm','1cm','location','northeast')
end

xlabel('\textbf{Age of oceanic lithosphere (Ma)}','Interpreter','latex','FontSize',22)
ylabel('$\mathbf{\Delta t^{\ast}_S}$','Interpreter','latex','FontSize',22)
set(gca,'ylim',[-1 4],'xlim',[0 12],'fontsize',16,'box','on','linewidth',2)
set(gca,'ytick',[-1:1:4])


if ifsave
if ifpreds, prstr = '_wpreds'; else prstr = ''; end
save2pdf(31,['TRENDS_OBS_P_dtstar_vs_age',prstr],'figs')
save2pdf(32,['TRENDS_OBS_S_dtstar_vs_age',prstr],'figs')
end

%% ==================================================================== %% 
%% =============================    dT    ============================= %% 
%% ==================================================================== %% 

% % LOAD RESULTS
% load([resdir,'all_dT_P_Z.mat']); Pall_dT = all_dT;
% load([resdir,'all_dT_S_T.mat']); Sall_dT = all_dT;
% 
% [ Psta_dT] = lsq_sta_evt( Pall_dT,0.01);
% Pstd = nanvariance(Pall_dT,false,2).^(1/2);
% PNobs = sum(~isnan(Pall_dT),2);
% 
% [ Ssta_dT] = lsq_sta_evt( Sall_dT,0.01);
% Sstd = nanvariance(Sall_dT,false,2).^(1/2);
% SNobs = sum(~isnan(Sall_dT),2);

%% Work out good stations, and plot them
%P
Pgd  = dTP.dT~=0 & abs(dTP.dT)<5 & dTP.Nobs>=5;
gdjdf  = Pgd & dTP.isob & dTP.ijdf;
gdgor  = Pgd & dTP.isob & ~dTP.ijdf;
dcshift = 0.3;
figure(33),clf, set(gcf,'pos',[400 500 800 360]),hold on
if ifpreds
for ig = 1:length(gs)
hp(ig) = plot(ages,dtt_P(:,ig)-min(dtt_P(:,ig)),'k','LineWidth',2,'LineStyle',gsls{ig});
end
end
scatter(abs(dTP.Xrdg(gdgor)/km_per_ma),dTP.dT(gdgor)+dcshift,20*dTP.Nobs(gdgor),...
    'k','MarkerFaceColor','b','LineWidth',1.5)
scatter(abs(dTP.Xrdg(gdjdf)/km_per_ma),dTP.dT(gdjdf)+dcshift,20*dTP.Nobs(gdjdf),...
    'k','MarkerFaceColor','r','LineWidth',1.5)
if ifpreds
legend(hp,'0.1mm','1mm','1cm','location','northeast')
end

xlabel('\textbf{Age of oceanic lithosphere (Ma)}','Interpreter','latex','FontSize',22)
ylabel('$\mathbf{\delta T_P}$','Interpreter','latex','FontSize',22)
set(gca,'ylim',[-0.1 1.5],'xlim',[0 12],'fontsize',16,'box','on','linewidth',2)


%S
Sgd  = dTS.dT~=0 & abs(dTS.dT)<5 & dTS.Nobs>=5;
gdjdf  = Sgd & dTS.isob & dTS.ijdf;
gdgor  = Sgd & dTS.isob & ~dTS.ijdf;
dcshift = 0.5;
figure(34),clf, set(gcf,'pos',[400 000 800 360]),hold on
if ifpreds
for ig = 1:length(gs)
hs(ig) = plot(ages,dtt_S(:,ig)-min(dtt_S(:,ig)),'k','LineWidth',2,'LineStyle',gsls{ig});
end
end
scatter(abs(dTS.Xrdg(gdgor)/km_per_ma),dTS.dT(gdgor)+dcshift,20*dTS.Nobs(gdgor),...
    'k','MarkerFaceColor','b','LineWidth',1.5)
scatter(abs(dTS.Xrdg(gdjdf)/km_per_ma),dTS.dT(gdjdf)+dcshift,20*dTS.Nobs(gdjdf),...
    'k','MarkerFaceColor','r','LineWidth',1.5)
if ifpreds
legend(hs,'0.1mm','1mm','1cm','location','northeast')
end

xlabel('\textbf{Age of oceanic lithosphere (Ma)}','Interpreter','latex','FontSize',22)
ylabel('$\mathbf{\delta T_S}$','Interpreter','latex','FontSize',22)
set(gca,'ylim',[-1 4],'xlim',[0 12],'fontsize',16,'box','on','linewidth',2)
set(gca,'ytick',[-1:1:4])

if ifsave
if ifpreds, prstr = '_wpreds'; else prstr = ''; end
save2pdf(33,['TRENDS_OBS_P_dT_vs_age',prstr],'figs')
save2pdf(34,['TRENDS_OBS_S_dT_vs_age',prstr],'figs')
end