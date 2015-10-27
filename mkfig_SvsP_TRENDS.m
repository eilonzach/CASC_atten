clear all
close all
cd ~/Documents/MATLAB/CASC_atten/
addpath('matguts')

%% PARMS

%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascattendb';
% DATA DIRECTORY (top level)
datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash
% RESULTS DIRECTORY
resdir = '~/Documents/MATLAB/CASC_atten/results/'; % needs final slash


%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%
% LOAD RESULTS
load([resdir,'all_dT_P_Z.mat']); dTp = all_dT;
load([resdir,'all_dT_S_T.mat']); dTs = all_dT;
load([resdir,'all_dtstar_P_Z.mat']); dtstp = all_dtstar;
load([resdir,'all_dtstar_S_T.mat']); dtsts = all_dtstar;
[ nstas,stas,slats,slons,selevs,ondate,offdate,staname ] = db_stadata( dbdir,dbnam );

%% ========================  STATION AVERAGES  ========================= %%
[ stav_dTp,~ ] = lsq_sta_evt( dTp,0.01 );
[ stav_dTs,~ ] = lsq_sta_evt( dTs,0.01 );
[ stav_dtstp,~ ] = lsq_sta_evt( dtstp,0.01 );
[ stav_dtsts,~ ] = lsq_sta_evt( dtsts,0.01 );

stav_ydT = stav_dTp~=0 & stav_dTs~=0;
stav_ydts = stav_dTp~=0 & stav_dTs~=0;

% stav diff-tt line fit
tts_stav = stav_dTs(stav_ydT) - mean( stav_dTs(stav_ydT) );
ttp_stav = stav_dTp(stav_ydT) - mean( stav_dTp(stav_ydT) );
[ m_tt_stav,~,mstd,bstd ] = fit_LSqEr( tts_stav,ttp_stav,1,1,1:0.01:5,[],10,1);

% stav diff-tstar line fit
tss_stav = stav_dtsts(stav_ydts) - mean( stav_dtsts(stav_ydts) );
tsp_stav = stav_dtstp(stav_ydts) - mean( stav_dtstp(stav_ydts) );
% [ m_ts,b,mstd,bstd ] = fit_LSqEr( tsp,tss,1,1,1:0.01:5,[],4,1);
figure(31);scatter(slons(stav_ydts),slats(stav_ydts),100,tss_stav,'filled')
figure(31);scatter(slons(stav_ydtp),slats(stav_ydtp),100,tsp_stav,'filled')
colormap(parula)
return

%% PLOT STAV DIFFERENTIAL TRAVEL TIME

figure(11), clf, set(gcf,'position',[200 200 480 610])
hold on

scatter(ttp_stav,tts_stav,40,'b','o','LineWidth',1.5)

axis([-3.7 3.7 -5 5]); ax = axis;

plot(ax(1:2),m_tt_stav*ax(1:2),'r--','LineWidth',3)

text(-3.2,4.6,sprintf('\\textbf{Slope = %.2f}',m_tt_stav),'FontSize',17,'interpreter','latex')
text(-3,4.15,sprintf('\\textbf{$\\Rightarrow\\, \\delta \\ln{V_S}/\\delta \\ln{V_P} = %.2f$}',m_tt_stav/1.81),'FontSize',17,'interpreter','latex')

set(gca,'fontsize',12,'color','none')
% title('Ratio of $P$ and $S$ travel times residuals','FontSize',18,'interpreter','latex')
xlabel('Sta-Av $P$ residual time (s)','FontSize',16,'interpreter','latex')
ylabel('Sta-Av $S$ residual time (s)','FontSize',16,'interpreter','latex')
return
save2pdf(11,'dTs_vs_dTp_stav','figs')

%% PLOT STAV DIFFERENTIAL T-STAR
figure(12), clf, set(gcf,'position',[300 200 480 610])
hold on

scatter(tsp_stav,tss_stav,40,'b','o','LineWidth',1.5)

axis([-3.7 3.7 -5 5]); ax = axis;

plot(ax(1:2),m_ts*ax(1:2),'r--','LineWidth',3)

text(-3.2,4.6,sprintf('\\textbf{Slope = %.2f}',m_ts),'FontSize',17,'interpreter','latex')
text(-3,4.15,sprintf('\\textbf{$\\Rightarrow\\, Q_P/Q_S = %.2f$}',m_ts/1.81),'FontSize',17,'interpreter','latex')

set(gca,'fontsize',12,'color','none')
% title('Ratio of $P$ and $S$ travel times residuals','FontSize',18,'interpreter','latex')
xlabel('$P$-wave $\Delta t^*$ (s)','FontSize',16,'interpreter','latex')
ylabel('$S$-wave $\Delta t^*$ (s)','FontSize',16,'interpreter','latex')

save2pdf(12,'dtstars_vs_dtstarp','figs')




%% ======================  INDIVIDUAL MEASUREMENTS  ====================== %%
ydT = ~isnan(dTp) & ~isnan(dTs);
ydts = ~isnan(dtstp) & ~isnan(dtsts);
yS = ~isnan(dtsts) & ~isnan(dTs);
yP = ~isnan(dtstp) & ~isnan(dTp);

% diff-tt line fit
tts = dTs(ydT) - mean( dTs(ydT) );
ttp = dTp(ydT) - mean( dTp(ydT) );
[ m_tt,b,mstd,bstd ] = fit_LSqEr( ttp,tts,1,1,1:0.01:5,[],4,1);

% diff-tstar line fit
tss = dtsts(ydts) - mean( dtsts(ydts) );
tsp = dtstp(ydts) - mean( dtstp(ydts) );
[ m_ts,b,mstd,bstd ] = fit_LSqEr( tsp,tss,1,1,1:0.01:5,[],4,1);



%% PLOT DIFFERENTIAL TRAVEL TIME

figure(1), clf, set(gcf,'position',[200 200 480 610])
hold on

scatter(ttp,tts,40,'b','o','LineWidth',1.5)

axis([-3.7 3.7 -5 5]); ax = axis;

plot(ax(1:2),m_tt*ax(1:2),'r--','LineWidth',3)

text(-3.2,4.6,sprintf('\\textbf{Slope = %.2f}',m_tt),'FontSize',17,'interpreter','latex')
text(-3,4.15,sprintf('\\textbf{$\\Rightarrow\\, \\delta \\ln{V_S}/\\delta \\ln{V_P} = %.2f$}',m_tt/1.81),'FontSize',17,'interpreter','latex')

set(gca,'fontsize',12,'color','none')
% title('Ratio of $P$ and $S$ travel times residuals','FontSize',18,'interpreter','latex')
xlabel('$P$ residual time (s)','FontSize',16,'interpreter','latex')
ylabel('$S$ residual time (s)','FontSize',16,'interpreter','latex')

save2pdf(1,'dTs_vs_dTp','figs')


%% PLOT DIFFERENTIAL T-STAR
figure(2), clf, set(gcf,'position',[300 200 480 610])
hold on

scatter(tsp,tss,40,'b','o','LineWidth',1.5)

axis([-3.7 3.7 -5 5]); ax = axis;

plot(ax(1:2),m_ts*ax(1:2),'r--','LineWidth',3)

text(-3.2,4.6,sprintf('\\textbf{Slope = %.2f}',m_ts),'FontSize',17,'interpreter','latex')
text(-3,4.15,sprintf('\\textbf{$\\Rightarrow\\, Q_P/Q_S = %.2f$}',m_ts/1.81),'FontSize',17,'interpreter','latex')

set(gca,'fontsize',12,'color','none')
% title('Ratio of $P$ and $S$ travel times residuals','FontSize',18,'interpreter','latex')
xlabel('$P$-wave $\Delta t^*$ (s)','FontSize',16,'interpreter','latex')
ylabel('$S$-wave $\Delta t^*$ (s)','FontSize',16,'interpreter','latex')

save2pdf(2,'dtstars_vs_dtstarp','figs')


%% PLOT S tt vs. tstar
figure(3)
plot(dTs(yS),dtsts(yS),'o')
title('Ratio of differential travel time to t-star for S','FontSize',18,'interpreter','latex')
xlabel('$\delta T_S$','FontSize',15,'interpreter','latex')
ylabel('$\Delta t^*_S$','FontSize',15,'interpreter','latex')

%% PLOT P tt vs. tstar
figure(4)
plot(dTp(yP),dtstp(yP),'o')
title('Ratio of differential travel time to t-star for P','FontSize',18,'interpreter','latex')
xlabel('$\delta T_P$','FontSize',15,'interpreter','latex')
ylabel('$\Delta t^*_P$','FontSize',15,'interpreter','latex')


