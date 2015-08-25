addpath('synthetic_tstar')
%% My simple synthetics
redo = 0;
pred_tstar

%% Goes 2012 synthetics
Goes12_synth

%%  ================= Master plot with all =================
figure(13), clf, set(gcf,'position',[0 0 1500 800])
colormap(flipud(jet))
% plot Q
subplot(4,5,1:2)
contourf(ages,Z,Q_s(:,:,2),[0:5:300]); shading flat
set(gca,'Ydir','reverse','FontSize',10,'XTickLabel',[])
set(gca,'position',[0.1300    0.7673    0.2866    0.1577])
set(get(gca,'Ylabel'),'string','Depth (km)','Interpreter','Latex','FontSize',12,'Position',[-2.2 102.5  1.0001])
set(get(gca,'Title'),'string','S-waves','Interpreter','Latex','FontSize',14)
text(0.4,0.85*200,'a)','FontSize',15), box on

subplot(4,5,3:4)
contourf(ages,Z,Q_p(:,:,2),[0:5:300]); shading flat
set(gca,'Ydir','reverse','FontSize',10,'XTickLabel',[],'YTickLabel',[])
set(gca,'position',[0.4556    0.7673    0.2866    0.1577])
set(get(gca,'Title'),'string','P-waves','Interpreter','Latex','FontSize',14)
text(0.4,0.85*200,'b)','FontSize',15), box on

subplot(4,5,5)
hc = colorbar; caxis([0 300])
set(gca,'visible','off')
set(gca,'position',[0.7813    0.7673    0.1    0.1577])
set(hc,'position',[0.77 0.7665 0.018 0.1583],'FontSize',10)
set(get(hc,'ylabel'),'string','Quality factor','Interpreter','Latex','FontSize',12,'position',[9 130 1])

% plot dV
subplot(4,5,6:7)
contourf(ages,Z,100*(V_s(:,:,2)-V_s_av(:,:,2))./V_s_av(:,:,2),[-10,-5:0.2501:5]); shading flat
set(gca,'Ydir','reverse','FontSize',10,'XTickLabel',[])
set(gca,'position',[0.1300    0.5482    0.2866    0.1577])
set(get(gca,'Ylabel'),'string','Depth (km)','Interpreter','Latex','FontSize',12,'Position',[-2.2 102.5  1.0001])
text(0.4,0.85*200,'c)','FontSize',15), box on

subplot(4,5,8:9)
contourf(ages,Z,100*(V_p(:,:,2)-V_p_av(:,:,2))./V_p_av(:,:,2),[-10,-5:0.2501:5]); shading flat
set(gca,'Ydir','reverse','FontSize',10,'XTickLabel',[],'YTickLabel',[])
set(gca,'position',[0.4556    0.5482    0.2866    0.1577])
text(0.4,0.85*200,'d)','FontSize',15), box on

subplot(4,5,10)
hc = colorbar; caxis([-10 5])
set(gca,'visible','off')
set(gca,'position',[0.7813    0.5482    0.1    0.1577])
set(hc,'position',[0.77 0.5475 0.018 0.1583])
set(get(hc,'ylabel'),'string','$\delta V_{S,P}~$ (\%)','Interpreter','Latex','FontSize',13,'position',[9 -3 1])

% plot dtstar_S,_P
subplot(4,5,11:12), hold on
for ig = 1:length(gs)
plot(ages,dtstar_S(:,ig),'r','LineWidth',2,'LineStyle',gsls{ig})
end
grid on, box on
set(gca,'FontSize',10,'XtickLabel',[],'Ylim',[-0.7 0],'XLim',[0 max(ages)])
set(gca,'position',[0.1300    0.3291    0.2866    0.1577])
ylabel('$\Delta t^{\star}_S$','Interpreter','latex','FontSize',14,'Position',[-2.2 -0.3587 1.0001])
plot(10,dtstar_dhyd,'or','MarkerSize',7,'LineWidth',2)
plot(10,dtstar_undhyd,'sr','MarkerSize',7,'LineWidth',2)
text(0.4,0.85*-0.7,'e)','FontSize',15)

subplot(4,5,13:14), hold on
for ig = 1:length(gs)
plot(ages,dtstar_P(:,ig),'r','LineWidth',2,'LineStyle',gsls{ig})
end
grid on, box on
set(gca,'FontSize',10,'XtickLabel',[],'YAxisLocation','right','Ylim',[-0.2 0],'XLim',[0 max(ages)])
set(gca,'position',[0.4556    0.3291    0.2866    0.1577])
ylabel('$\Delta t^{\star}_P$','Interpreter','latex','FontSize',14,'Position',[17.2 -0.1025  1.0001])
text(0.4,0.85*-0.2,'f)','FontSize',15)

% plot dtt_S,_P
subplot(4,5,16:17), hold on
for ig = 1:length(gs)
plot(ages,dtt_S(:,ig),'b','LineWidth',2,'LineStyle',gsls{ig})
end
grid on, box on
set(gca,'FontSize',10,'XLim',[0 max(ages)],'YLim',[-2 0])
set(get(gca,'xlabel'),'string','Age of oceanic plate (Ma)','Interpreter','latex','FontSize',12,'position',[7. -2.4 1.00011])
set(gca,'position',[0.1300    0.11    0.2866    0.1577])
ylabel('$\Delta t_S$','Interpreter','latex','FontSize',14,'Position',[-2.2 -1.02  1.0001])
plot(10,dT_dhyd,'ob','MarkerSize',7,'LineWidth',2)
plot(10,dT_undhyd,'sb','MarkerSize',7,'LineWidth',2)
text(0.4,0.85*-2,'g)','FontSize',15)

subplot(4,5,18:19), hold on
for ig = 1:length(gs)
plot(ages,dtt_P(:,ig),'b','LineWidth',2,'LineStyle',gsls{ig})
end
grid on, box on
set(gca,'FontSize',10,'YAxisLocation','right','XLim',[0 max(ages)],'YLim',[-1 0])
set(get(gca,'xlabel'),'string','Age of oceanic plate (Ma)','Interpreter','latex','FontSize',12,'position',[7. -1.2 1.00011])
set(gca,'position',[0.4556    0.11    0.2866    0.1577])
ylabel('$\Delta t_P$','Interpreter','latex','FontSize',14,'Position',[17.2 -0.51  1.0001])
text(0.4,0.85*-1,'h)','FontSize',15)


save2pdf(13,'synthetic_dt_dtstar','figs')