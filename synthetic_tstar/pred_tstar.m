% script to predict t* at a number of ages etc.
% ASSUMING CONSTANT CUSTOM GRAINSIZES

if ~exist('redo','var'), redo = 1; end
if redo

ages = [0:1:15]'; % age in My
Z = [00:5:200]';
gs = [1e-3 1e-2, 1e-1];
gsls = {':','--','-'};
frq = 1/10;

%% Loop and calculate
tstar_S = zeros(length(ages),length(gs));
tt_S = zeros(length(ages),length(gs));
tstar_P = zeros(length(ages),length(gs));
tt_P = zeros(length(ages),length(gs));

Q_s = zeros(length(Z),length(ages),length(gs));
Q_P = zeros(length(Z),length(ages),length(gs));
V_s = zeros(length(Z),length(ages),length(gs));
V_P = zeros(length(Z),length(ages),length(gs));
for ig = 1:length(gs)
for ia = 1:length(ages)

[ Qs1,Qp,Vs1,Vp ] = QV_at_z( ages(ia),Z,gs(ig),frq);
% [ Qs1,Qp,Vs1,Vp ] = QV_at_z( ages(ia),Z,gs(ig),frq);
[ Qs2,Qp2,Vs2,Vp2,gs_i ] = QV_at_z_WET( ages(ia),Z,1e2,frq);
Qs = Qs1;
Vs = Vs1;
if ages(ia)>14
figure(1), clf, 
subplot(121), hold on
plot(Vs1,Z);
plot(Vs2,Z,'r');
xlim([0,400]),set(gca,'ydir','reverse')
subplot(122), hold on
plot(gs(ig)*ones(size(Z)),Z);
plot(gs_i,Z,'r');
set(gca,'ydir','reverse','xscale','log')
xlim([1e-6,1])

end
% Qs(Qs>1000) = 1000;
Q_s(:,ia,ig) = Qs;
Q_p(:,ia,ig) = Qp;
V_s(:,ia,ig) = Vs;
V_p(:,ia,ig) = Vp;

for ii = 1:length(Z)-1
    dz = diff(Z(ii:ii+1))*1000;

    qsav = 1/(mean(Qs(ii:ii+1)));
    vsav = mean(Vs(ii:ii+1));
    tstar_S(ia,ig) = tstar_S(ia,ig) + qsav*dz/vsav;
    tt_S(ia,ig) = tt_S(ia,ig) + dz/vsav;
    
    qpav = 1/(mean(Qp(ii:ii+1)));
    vpav = mean(Vp(ii:ii+1));
    tstar_P(ia,ig) = tstar_P(ia,ig) + qpav*dz/vpav;
    tt_P(ia,ig) = tt_P(ia,ig) + dz/vpav;
    
    
    
end % tstar integration loop

end % loop on ages
end % loop on grain size

for ig = 1:length(gs)
    V_s_av(:,:,ig) = mean(V_s(:,:,ig),2)*ones(1,length(ages));
    V_p_av(:,:,ig) = mean(V_p(:,:,ig),2)*ones(1,length(ages));
end

% return

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
for ig = 1:length(gs)
plot(ages,dtstar_S(:,ig),'r',ages,dtstar_P(:,ig),'b','LineWidth',2,'LineStyle',gsls{ig})
end
xlabel('Age of oceanic lithosphere (Ma)','FontSize',15)
ylabel('$\Delta t^{\star}$','Interpreter','latex','FontSize',15)

subplot(212), hold on
for ig = 1:length(gs)
plot(ages,dtt_S(:,ig),'r',ages,dtt_P(:,ig),'b','LineWidth',2,'LineStyle',gsls{ig})
end
xlabel('Age of oceanic lithosphere (Ma)','FontSize',15)
ylabel('$\Delta t (s)$','Interpreter','latex','FontSize',15)

end % if redo




