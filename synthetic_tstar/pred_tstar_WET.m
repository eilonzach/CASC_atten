% script to predict t* at a number of ages etc.
% WITH GRAINSIZES CONTROLLED BY PALEOWATTMETER AND SOME WATER

if ~exist('redo','var'), redo = 1; end
if redo
% clear all

ages = [0:1:15]'; % age in My
Z = [00:5:200]';
ww = [50 200 1000]; %in wtPPM
wwls = {':','--','-'};
frq = 1/10;

%% Loop and calculate
tstar_S = zeros(length(ages),length(ww));
tt_S = zeros(length(ages),length(ww));
tstar_P = zeros(length(ages),length(ww));
tt_P = zeros(length(ages),length(ww));

Q_s = zeros(length(Z),length(ages),length(ww));
Q_P = zeros(length(Z),length(ages),length(ww));
V_s = zeros(length(Z),length(ages),length(ww));
V_P = zeros(length(Z),length(ages),length(ww));
for iw = 1:length(ww)
for ia = 1:length(ages)

[ Qs,Qp,Vs,Vp ] = QV_at_z_WET( ages(ia),Z,ww(iw),frq);
% Qs(Qs>1000) = 1000;
Q_s(:,ia,iw) = Qs;
Q_p(:,ia,iw) = Qp;
V_s(:,ia,iw) = Vs;
V_p(:,ia,iw) = Vp;

for ii = 1:length(Z)-1
    dz = diff(Z(ii:ii+1))*1000;

    qsav = 1/(mean(Qs(ii:ii+1)));
    vsav = mean(Vs(ii:ii+1));
    tstar_S(ia,iw) = tstar_S(ia,iw) + qsav*dz/vsav;
    tt_S(ia,iw) = tt_S(ia,iw) + dz/vsav;
    
    qpav = 1/(mean(Qp(ii:ii+1)));
    vpav = mean(Vp(ii:ii+1));
    tstar_P(ia,iw) = tstar_P(ia,iw) + qpav*dz/vpav;
    tt_P(ia,iw) = tt_P(ia,iw) + dz/vpav;
    
    
    
end % tstar integration loop

end % loop on ages
end % loop on grain size

for iw = 1:length(ww)
    V_s_av(:,:,iw) = mean(V_s(:,:,iw),2)*ones(1,length(ages));
    V_p_av(:,:,iw) = mean(V_p(:,:,iw),2)*ones(1,length(ages));
end



figure(11), clf
% plot Q
subplot(2,2,1)
contourf(ages,Z,Q_s(:,:,2),[0:10:300])
set(gca,'Ydir','reverse'); shading flat
subplot(2,2,2)
contourf(ages,Z,Q_p(:,:,2),[0:10:300])
set(gca,'Ydir','reverse'); shading flat
% % plot V
% subplot(2,2,3)
% contourf(ages,Z,V_s(:,:,2),[4.2:0.05:4.6]*1e3)
% set(gca,'Ydir','reverse'); shading flat
% subplot(2,2,4)
% contourf(ages,Z,V_p(:,:,2),[7.4:0.05:8.5]*1e3)
% set(gca,'Ydir','reverse'); shading flat
% plot dV
subplot(2,2,3)
contourf(ages,Z,100*(V_s(:,:,2)-V_s_av(:,:,2))./V_s_av(:,:,2),[-10,-5:1.001:5.015])
set(gca,'Ydir','reverse'); shading flat
subplot(2,2,4)
contourf(ages,Z,100*(V_p(:,:,2)-V_p_av(:,:,2))./V_p_av(:,:,2),[-10,-5:1.001:5.015])
set(gca,'Ydir','reverse'); shading flat

%% differential tt/tstar
dtstar_S = tstar_S - ones(length(ages),1)*tstar_S(1,:);
dtstar_P = tstar_P - ones(length(ages),1)*tstar_P(1,:);
dtt_S = tt_S - ones(length(ages),1)*tt_S(1,:);
dtt_P = tt_P - ones(length(ages),1)*tt_P(1,:);

%% plot tt and tstar
figure(12), clf
subplot(211), hold on
for iw = 1:length(ww)
plot(ages,dtstar_S(:,iw),'r',ages,dtstar_P(:,iw),'b','LineWidth',2,'LineStyle',wwls{iw})
end
xlabel('Age of oceanic lithosphere (Ma)','FontSize',15)
ylabel('$\Delta t^{\star}$','Interpreter','latex','FontSize',15)

subplot(212), hold on
for iw = 1:length(ww)
plot(ages,dtt_S(:,iw),'r',ages,dtt_P(:,iw),'b','LineWidth',2,'LineStyle',wwls{iw})
end
xlabel('Age of oceanic lithosphere (Ma)','FontSize',15)
ylabel('$\Delta t (s)$','Interpreter','latex','FontSize',15)

end % if redo




