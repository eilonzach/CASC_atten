function [misfit_amp,misfit_phi,E_a,E_p] = plot_AmpPhi_fit( Amp,Phi,fmids,wts,dtstar,dT,A0, alp,fign )
% [misfit_amp,misfit_phi] = plot_AmpPhi_fit( Amp,Phi,fmids,[wts=1],dtstar,dT,A0, [alpha=0], [fign=33] )
% 
% Quick function to plot amplitude and phase spectra and the fits to them,
% along with the misfit.

if nargin < 8 || isempty(alp)
    alp = 0;
end
if nargin < 9 || isempty(fign)
    fign = 33;
end
if isempty(wts)
    wts = ones(size(Amp));
end

Amp = Amp(:);
Phi = Phi(:);
fmids = fmids(:);
ws = 2*pi*fmids;

ffa = fmids;
ffp = fmids;
xlaba = 'freq (Hz)';
xlabq = 'freq (Hz)';

if alp == 0;
    lnApred   = log(A0) - pi*fmids*dtstar;
    phipred = -(1/pi)*log(fmids)*dtstar + dtstar + dT;
%     ffa = fmids;
%     ffp = fmids;
%     xlaba = 'freq (Hz)';
%     xlabq = 'freq (Hz)';
else
    lnApred = log(A0) - 0.5*ws.^(1-alp)*dtstar;
    phipred = 0.5*ws.^(-alp).*(2*pi).^alp * cot(alp*pi/2)*dtstar + dT;
%     ffa = fmids.^(1-alp);
%     ffp = fmids.^(-alp);
%     xlaba = 'freq$^{(1-\alpha)}$';
%     xlabq = 'freq$^{(-\alpha)}$';
end

E_a = Amp-exp(lnApred);
E_p = Phi-phipred;

misfit_amp = E_a'*diag(wts)*E_a;
misfit_phi = E_p'*diag(wts)*E_p;


%% START THE FANCY PLOTS
figure(fign), clf, set(gcf,'pos',[600 600 600,800])


% ================ PLOT THE AMPLITUDE SPECTRA ================
subplot(2,1,1), hold on

scatter(ffa,log(Amp),150*wts,'ok','MarkerFaceColor','r','linewidth',1.5)
plot(ffa,lnApred,'g','Linewidth',2.5)

xlabel(xlaba,'FontSize',22,'interpreter','latex')
ylabel('$\ln\,(R_{12})$','FontSize',22,'interpreter','latex')
title(sprintf('$\\Delta t^*$ = %.2f \\,\\,\\, $\\alpha$ = %.2f',dtstar,alp),'FontSize',22,'interpreter','latex')
set(gca,'fontsize',15,'Xscale','linear','ylim',[-3 1],'linewidth',2,'box','on')

ax = axis; xl = ax(1) + 0.05*(ax(2)-ax(1)); yl = ax(4) - 0.1*(ax(4)-ax(3));
text(xl,yl,['Weighted misfit = ',num2str(misfit_amp)],'interpreter','latex','fontsize',20,'horizontalalignment','left','verticalalignment','bottom')


% ================ PLOT THE PHASE SPECTRA ================
subplot(2,1,2), hold on
scatter(ffp,Phi,150*wts,'ok','MarkerFaceColor','r','linewidth',1.5)
plot(ffp,phipred,'g','Linewidth',2.5)

% ha1 = arrow([eqar(ind1).fcross,1],[eqar(ind1).fcross,0],5,90,2,0.1,'FaceColor','m'); % fcross for sta 1
% ha2 = arrow([eqar(ind2).fcross,1],[eqar(ind2).fcross,0],5,90,2,0.1,'FaceColor','g'); % fcross for sta 2
% plot(fmax*[1 1],[-1 1],'--b')  

xlabel(xlabq,'FontSize',22,'interpreter','latex')
ylabel('$\Delta \phi_{12}$ (s)','FontSize',22,'interpreter','latex')
set(gca,'fontsize',15,'Xscale','linear','linewidth',2,'box','on')

ax = axis; xl = ax(1) + 0.05*(ax(2)-ax(1)); yl = ax(4) - 0.1*(ax(4)-ax(3));
text(xl,yl,['Weighted misfit = ',num2str(misfit_phi)],'interpreter','latex','fontsize',20,'horizontalalignment','left','verticalalignment','bottom')

return
% ================ SAVE FIGURE ================
ofile = sprintf('synth_comb_Q1-%.0f_Q2-%.0f',Q0_1,Q0_2);
save2pdf(33,ofile,'figs');


end

