clear all
L = [5:5:200];
Q = logspace(0,2,40);   
V = 4.1;                    

Nq = length(Q);
Nl = length(L);

tstlim = [0:0.5:2.5];
Qlim = [1 70];

tstar = (L'*ones(1,Nq))./V./(ones(Nl,1)*Q);

figure(9); clf, hold on
contourf((ones(Nl,1)*Q),(L'*ones(1,Nq)),tstar,[tstlim(1):0.05:tstlim(end)])
shading flat
contour((ones(Nl,1)*Q),(L'*ones(1,Nq)),tstar,tstlim,'k','LineWidth',1.5)

text(max(Qlim)-2,min(L)+15,sprintf('V_S = %.1f km/s',V),...
    'color','w','Fontsize',18,'FontWeight','bold','HorizontalAlignment','right')

hc = colorbar;
caxis([tstlim(1),tstlim(end)])

set(gca,'xlim',Qlim,'box','on','Fontsize',14,'LineWidth',2)
xlabel('$\mathbf{Q}$','Fontsize',20,'interpreter','latex')
ylabel('$\mathbf{L \,\,(km)}$','Fontsize',20,'interpreter','latex')

% put on box over it
ax1 = gca;
ax2 = axes;
set(ax2,'position',get(ax1,'position'),...
    'xtick',get(ax1,'xtick'),'xticklabel',[],...
    'ytick',get(ax1,'ytick'),'yticklabel',[],...
    'box','on','linewidth',2,...
    'color','none')


save2pdf(9,'tstar_Q_L_function','/Users/zeilon/Documents/MATLAB/CASC_atten/figs')