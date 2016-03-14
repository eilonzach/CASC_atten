% script to make some snthetic dT and dtstar predictions based on the
% models of Goes et al., JGR, 2012

%% 1315D case, with dehyd,depl,melt (Figure 5d)
[Vs_0,Zv_0]   = textread('Goes12_fig5d_0Ma_Vs.txt','%f %f','headerlines',1);
[Vs_10,Zv_10] = textread('Goes12_fig5d_10Ma_Vs.txt','%f %f','headerlines',1);
[Qs_0,Zq_0]   = textread('Goes12_fig6a_0Ma_dhyd_Qs.txt','%f %f','headerlines',1);
[Qs_10,Zq_10] = textread('Goes12_fig6a_10Ma_dhyd_Qs.txt','%f %f','headerlines',1);


Z_g = [1:1:180]';
Vs_g = zeros(length(Z_g),2);
Qs_g = zeros(length(Z_g),2);
Vs_g(:,1) = linterp(Zv_0,Vs_0,Z_g);
Vs_g(:,2) = linterp(Zv_10,Vs_10,Z_g);
Qs_g(:,1) = linterp(Zq_0,Qs_0,Z_g);
Qs_g(:,2) = linterp(Zq_10,Qs_10,Z_g);

dz_g = [diff(Z_g),diff(Z_g)];
Vsav_g = 0.5*(Vs_g(1:end-1,:)+Vs_g(2:end,:));
Qsav_g = 0.5*(Qs_g(1:end-1,:)+Qs_g(2:end,:));

tstar_g = sum(dz_g./Vsav_g./Qsav_g);
tt_g = sum(dz_g./Vsav_g);

dtstar_dhyd = diff(tstar_g);
dT_dhyd = diff(tt_g);

fprintf('dt*_s between 0 and 10 Ma for Goes12 dhyd model is %.2f\n',dtstar_dhyd)
fprintf('dT_s between 0 and 10 Ma for Goes12 dhyd model is %.2f\n',dT_dhyd)

%% 1315D case, without dehyd,depl,melt (Figure 5a)
[Vs_0,Zv_0]   = textread('Goes12_fig5a_0Ma_Vs.txt','%f %f','headerlines',1);
[Vs_10,Zv_10] = textread('Goes12_fig5a_10Ma_Vs.txt','%f %f','headerlines',1);
[Qs_0,Zq_0]   = textread('Goes12_fig6a_0Ma_undhyd_Qs.txt','%f %f','headerlines',1);
[Qs_10,Zq_10] = textread('Goes12_fig6a_10Ma_undhyd_Qs.txt','%f %f','headerlines',1);


Z_g = [1:1:180]';
Vs_g = zeros(length(Z_g),2);
Qs_g = zeros(length(Z_g),2);
Vs_g(:,1) = linterp(Zv_0,Vs_0,Z_g);
Vs_g(:,2) = linterp(Zv_10,Vs_10,Z_g);
Qs_g(:,1) = linterp(Zq_0,Qs_0,Z_g);
Qs_g(:,2) = linterp(Zq_10,Qs_10,Z_g);

dz_g = [diff(Z_g),diff(Z_g)];
Vsav_g = 0.5*(Vs_g(1:end-1,:)+Vs_g(2:end,:));
Qsav_g = 0.5*(Qs_g(1:end-1,:)+Qs_g(2:end,:));

tstar_g = sum(dz_g./Vsav_g./Qsav_g);
tt_g = sum(dz_g./Vsav_g);

dtstar_undhyd = diff(tstar_g);
dT_undhyd = diff(tt_g);

fprintf('dt*_s between 0 and 10 Ma for Goes12 undhyd model is %.2f\n',dtstar_undhyd)
fprintf('dT_s between 0 and 10 Ma for Goes12 undhyd model is %.2f\n',dT_undhyd)
