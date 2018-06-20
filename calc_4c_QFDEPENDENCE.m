% clear all
close all
cd /Users/zeilon/Documents/MATLAB/CASC_atten/

ifsave = false;

% orids = [124,140,142,170,196,238,263,269];
orids = 1:350;

phasecomp = 'ST';

ifobs = 0;

test_alphas = [0:0.05:0.9];
fmids = 1./logspace(log10(1),log10(20),30)';

amp2phiwt = 3;

norids = length(orids);
na = length(test_alphas);

%% re-do the fitting calculations
avr = zeros(na, norids);
aew = zeros(na, norids);
aps = zeros(norids,1);
for ie = 1:norids
    pairfile = sprintf('results_pairspecs/%.0f_pairspecs_%s.mat',orids(ie),phasecomp);
    if exist(pairfile,'file')~=2, continue, end
    load(pairfile)
    
    Amat = pairwise.As;
    phimat = pairwise.phis;
    wtmat = double(pairwise.inds).*pairwise.wts;
    
    if ifobs % pull out/use only OBS stations by setting wts to zero for pairs with any non-OBS
        nstas = length(sts);
        iob = zeros(nstas,1);
        hob = zeros(handshake(nstas),1);
        for is = 1:nstas, iob(is) = ~isempty(which_OBS(sts{is})); end
        k=0;
        for is1 = 1:nstas
        for is2 = is1+1:nstas
        k=k+1;
        wtmat(k,:) = wtmat(k,:)*double(iob(is1) & iob(is2));
        end
        end
    end


    [ delta_tstar_pref,delta_T_pref,A0_pref,alpha_pref,alpha_VR,alpha_Ew ] ...
        = calc_fdependent( Amat,phimat,fmids,test_alphas,wtmat,amp2phiwt,2,'',0);


    aps(ie,1) = alpha_pref;
    avr(:,ie) = alpha_VR;
    aew(:,ie) = alpha_Ew;
end


wts = zeros(norids,1);
for ie = 1:norids
    pairfile = sprintf('results_pairspecs/%.0f_pairspecs_ST.mat',orids(ie));
    if exist(pairfile,'file')~=2, continue, end
    load(pairfile)
    wts(ie) = sqrt(sum(sum(pairwise.wts)));
end

nsts = zeros(norids,1);
for ie = 1:norids
    pairfile = sprintf('results_pairspecs/%.0f_pairspecs_ST.mat',orids(ie));
    if exist(pairfile,'file')~=2, continue, end
    load(pairfile)
    nsts(ie) = quadratic_solve(1,-1,-2*length(pairwise.dtstar));
end
    


    
%% Calculate the normalised variance averaged over all events
aew2 = aew.^2; % square error to get variance
aewn = aew2./(ones(na,1)*min(aew2)); % normalise

nsts = nsts(~isnan(nanmean(aewn))); % chuck out nans
wtsn = wts(~isnan(nanmean(aewn))); % chuck out nans
aewn = aewn(:,~isnan(nanmean(aewn))); % chuck out nans

aewnALL = aewn*wtsn/sum(wtsn); % weighted average
aewnALL = aewnALL./min(aewnALL) % re-normalise
% F-test
fstatistic = fcdf(aewnALL,sum(wtsn),sum(wtsn));
[~,alp_sig2] = crossing(fstatistic,test_alphas,0.95); if isempty(alp_sig2), alp_sig2 = nan; end
[~,alp_sig1] = crossing(fstatistic,test_alphas,0.68); if isempty(alp_sig1), alp_sig1 = nan; end
sig2 = interp1(test_alphas,aewnALL,alp_sig2);
sig1 = interp1(test_alphas,aewnALL,alp_sig1);

%% plot
figure(76), clf, set(gcf,'pos',[100 100 570 400]), hold on
plot(test_alphas,aewnALL ,'-sk','Markerfacecolor','k','Linewidth',2);

xlabel('frequency exponent ($\alpha$)','interpreter','latex','FontSize',22)
ylabel('Relative misfit, ($\chi^2/\chi^2_{min}$)','interpreter','latex','FontSize',22)

set(gca,'FontSize',15,'box','on','ylim',.01*[99.8 102],'xlim',[0 max(test_alphas)],'FontName','Helvetica')

text(0.03,1.018,sprintf('%.0f events',length(nsts)),'Fontsize',18,'interpreter','latex')
text(0.03,1.0155,sprintf('%.1f stations (mean)',mean(nsts)),'Fontsize',18,'interpreter','latex')
text(0.03,1.013,sprintf('%.0f $\\Delta t^*$ measurements',sum(handshake(nsts))),'Fontsize',18,'interpreter','latex')

ax = axis;
plot([ax(1),alp_sig1,alp_sig1],[sig1 sig1,ax(3)],'--k','Linewidth',1)
text(0.03,sig1+0.001,'$1\sigma$','Fontsize',16,'interpreter','latex')
plot([ax(1),alp_sig2,alp_sig2],[sig2 sig2,ax(3)],'--k','Linewidth',1)
text(0.03,sig2+0.001,'$2\sigma$','Fontsize',16,'interpreter','latex')



if ifsave
    save2pdf(76,'alpha_tests_ALLwav','figs')
end
return

%% plot all individually
figure(77), clf;
h = plot(test_alphas,aewn ,'-o','MarkerFaceColor','auto');
xlabel('frequency exponent ($\alpha$)','interpreter','latex','FontSize',22)
ylabel('Relative misfit, ($\chi^2/\chi^2_{min}$)','interpreter','latex','FontSize',22)
set(gca,'FontSize',14,'box','on','ylim',.01*[99.8 102],'FontName','Helvetica')
for ih = 1:length(h)
    cl = colour_get(ih,length(h),1,parula);
    set(h(ih),'markerfacecolor',cl,'color',cl); 
end
hl = legend(num2str(orids'),'Location','northwest');
set(hl,'FontSize',12,'FontName','Helvetica');
if ifsave
    save2pdf(77,'alpha_tests','figs')
end
return

figure(79), clf, set(gcf,'pos',[10 10 500 800])
subplot(311)
plot(test_alphas,Eaw,'-or')
title([titlstr,' ',methstr],'FontSize',22,'interpreter','latex')
ylabel('Eaw','FontSize',18,'interpreter','latex')
set(gca,'XtickLabel',[],'xtick',test_alphas,'fontsize',14)

subplot(312)
plot(test_alphas,Epw,'-ob')
ylabel('Epw','FontSize',18,'interpreter','latex')
set(gca,'XtickLabel',[],'xtick',test_alphas,'fontsize',14)

subplot(313)
plot(test_alphas,alpha_misfits,'-ok')
ylabel('Wt misfit','FontSize',18,'interpreter','latex')
xlabel('alphas','FontSize',18,'interpreter','latex')
set(gca,'xtick',test_alphas,'fontsize',14)


%% sum using existing VR saved in structure
avr = zeros(19, norids);
for ie = 1:norids
    pairfile = sprintf('results_amisfits/%.0f_amisfits_ST.mat',orids(ie));
    if exist(pairfile,'file')~=2, continue, end
    load(pairfile)
    avr(:,ie) = alpha_VR;
end
figure(76), clf; hold on
plot([0:0.05:0.9],100*(2-mean(avr')) ,'-sk','Markerfacecolor','k','Linewidth',2);
xlabel('frequency exponent ($\alpha$)','interpreter','latex','FontSize',22)
ylabel('Relative misfit, ($\chi^2/\chi^2_{min}$)','interpreter','latex','FontSize',22)
set(gca,'FontSize',14,'box','on','ylim',.01*[99.8 107],'FontName','Helvetica')



