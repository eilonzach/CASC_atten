%% %% %% 
% clear all
saveopt = 1;
% ofile = '~/Documents/MATLAB/CASC_atten/TOMOGRAPHY/results/Ltest_dtstar_ST';
ofig = 'Lcurv_dtstar_ST';
figN = 6;
% load(ofile)


Nd = length(Ltest.damp);
Ns = length(Ltest.smooth);
cls_d = colour_get(1:Nd,Nd,1,flipud(autumn));
cls_s = colour_get(1:Ns,Ns+0.2,0.8,flipud(parula));
close all

%% Calc norms
norm = Ltest.norm + Ltest.norm_estt + Ltest.norm_sstt;
norm = Ltest.norm;

%% plot
figure(33), clf
set(gcf,'Position',[30 330 700 500])
ax1 = axes('position',[0.13 0.19 0.8 0.75]);
axes(ax1), hold on

for iss = 1:Ns
plot(Ltest.vr(:,iss),norm(:,iss),'-','MarkerSize',4,'Linewidth',4.,'color',cls_s(iss,:))
% scatter(Ltest.vr(:,iss),norm_all(:,iss),50,cls_d,'filled')
% scatter(Ltest.vr_tt(:,iss),norm_tt(:,iss),50,cls_d,'filled')
% scatter(Ltest.vr_dt(:,iss),norm_dT(:,iss),50,cls_d,'filled')  
for idd = 1:Nd
    plot(Ltest.vr(idd,iss),norm(idd,iss),'o','MarkerEdgeColor','none','MarkerSize',10,'MarkerFaceColor',cls_d(idd,:))
%     plot(Ltest.vr_tt(idd,iss),norm_tt(idd,iss),'o','MarkerEdgeColor','none','MarkerSize',5,'MarkerFaceColor',cls_d(idd,:))
%     plot(Ltest.vr_dt(idd,iss),norm_dT(idd,iss),'o','MarkerEdgeColor','none','MarkerSize',5,'MarkerFaceColor',cls_d(idd,:))
end
end

plot(Ltest.resid(Ltest.damp==3,Ltest.smooth==3),norm(Ltest.damp==3,Ltest.smooth==3),...
    'ok','MarkerSize',15,'MarkerFaceColor','k')

%% axis etc.
set(ax1,'Fontsize',16,'xdir','reverse','LineWidth',1.5,'XTick',[0:5:100],'YTick',[0:1:20])
xlabel('Variance reduction (\%)',  'FontSize',23,'Interpreter','Latex')
ylabel('Model norm (s)',          'FontSize',23,'Interpreter','Latex')
set(get(ax1,'Ylabel'),'Position',[100.8 3.5 1])
axis([85 100 0 7])

% text(101,-18,['\textbf{Figure ',num2str(figN),'}'],'Fontsize',25,'interpreter','latex')


%% legend + annotations
text(90.,99,'$\leftarrow$ less damped','interpreter','latex','fontsize',18,'rotation',-80)
text(62.8,19.4,'more damped $\rightarrow$ ','interpreter','latex','fontsize',18,'rotation',-8)
text(97,8,'$\gamma = 3$','interpreter','latex','fontsize',20)
text(96.7,14,'$\epsilon=3$','interpreter','latex','fontsize',20)
arrow([92,17],[87.2,34],'length',10)

%% key
keyloc = [90,5];
keysiz = [4,3];

kle = keyloc(1) - 0.5*keysiz(1); kri = keyloc(1) + 0.5*keysiz(1);
kbo = keyloc(2) - 0.5*keysiz(2); kto = keyloc(2) + 0.5*keysiz(2);
kw = keysiz(1); kh = keysiz(2);

hold on
patch([kle kle kri kri kle],[kbo kto kto kbo kbo],...
      'w','LineWidth',2)
text(keyloc(1),kto-0.1*kh,'smoothing ($\gamma$)',...
    'fontsize',15,'fontweight','bold','verticalalignment','top','horizontalalignment','center','interpreter','latex');

for iss = 1:Ns
plot(kle + [0.3*kw,0.85*kw],kbo + (iss/(Ns+3))*kh + [0 0],'LineWidth',2,'color',cls_s(iss,:))
text(kle + 0.2*kw,kbo + (iss/(Ns+3))*kh ,num2str(Ltest.smooth(iss)),'verticalalignment','middle','FontSize',12,'fontweight','bold')
end
% return
% scatter(kle + 0.3*kw,kbo + 0.61*kh,sqrt(1),  'k','filled','MarkerEdgeColor','k');
% scatter(kle + 0.3*kw,kbo + 0.435*kh,sqrt(10), 'k','filled','MarkerEdgeColor','k');
% scatter(kle + 0.3*kw,kbo + 0.18*kh,sqrt(50),'k','filled','MarkerEdgeColor','k');
% 
% text(kle + 0.62*kw,kbo + 0.61*kh,'1', 'fontsize',12,'fontweight','bold','verticalalignment','middle');
% text(kle + 0.62*kw,kbo + 0.435*kh,'10','fontsize',12,'fontweight','bold','verticalalignment','middle');
% text(kle + 0.62*kw,kbo + 0.18*kh,'50','fontsize',12,'fontweight','bold','verticalalignment','middle');




if saveopt
    save2pdf(33,ofig,'~/Documents/MATLAB/CASC_atten/TOMOGRAPHY/figs/');
end
