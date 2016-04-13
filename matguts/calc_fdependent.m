function [ delta_tstar_pref,delta_T_pref,A0_pref,alpha_pref,alpha_misfits ] ...
    = calc_fdependent( Amat,phimat,fmids,test_alphas,wtmat,amp2phiwt,opt,titlstr )
% [ delta_tstar_pref,delta_T_pref,A0_pref,alpha_pref ] 
%  = calc_fdependent( Amat,phimat,fmids,test_alphas,wtmat,[amp2phiwt=5],[opt=1],[titstr=''] )
%   Function to re-do all the calculations of dtstar and dT testing the
%   various frequency-dependent cases. 
% 
%  Opt==1 ==> all in one
%  Opt==2 ==> one-by-one
if nargin < 6 || isempty(amp2phiwt)
    amp2phiwt = 5;
end
if nargin < 7 || isempty(opt)
    opt = 2;
end
if nargin < 8 || isempty(titlstr)
    titlstr = '';
end

Na = length(test_alphas);
Eaw = zeros(Na,1);
Epw = zeros(Na,1);

if opt==1
 %% all in one method
    [ ~,~,~,~,~,dtstars_a,dTs_a,A0s_a ] ...
        = invert_allin1_Aphis_4_STA_dtdtstar_alpha(Amat,phimat,fmids,test_alphas,wtmat,amp2phiwt );
    methstr = 'one-by-one';
elseif opt == 2
   %% one-by-one method
    [ ~,~,~,~,~,dtstars_a,dTs_a,A0s_a ] ...
        = invert_1by1_Aphis_4_STA_dtdtstar_alpha(Amat,phimat,fmids,test_alphas,wtmat,amp2phiwt );
    methstr = 'all in one';
end

    

for ia = 1:length(test_alphas)
    % compute A and phi predictions for best fitting values at each alpha
    [ Amat_pred,phimat_pred ] = pred_Amat_phimat( dtstars_a(:,ia),dTs_a(:,ia),A0s_a(:,ia),fmids,test_alphas(ia) );
    Eaw(ia) = sum(sum((log(Amat_pred) - log(Amat)).^2.*wtmat));
    Epw(ia) = sum(sum((  phimat_pred -    phimat ).^2.*wtmat));
end
alpha_misfits = amp2phiwt*Eaw + Epw;


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

ia = mindex(alpha_misfits);

delta_tstar_pref = dtstars_a(:,ia);
delta_T_pref = dTs_a(:,ia);
A0_pref = A0s_a(:,ia);
alpha_pref = test_alphas(ia);


end

