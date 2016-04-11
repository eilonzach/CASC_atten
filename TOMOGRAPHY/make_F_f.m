function [ F,f,data ] = make_F_f( G,d,data,par )
% [ F,f,data ] = make_F_f( G,d,data,par )
% 
% create regularised, weighted matrices to invert, using the principle
% 
% [Wd * G]     [Wd*d]       Wd is the data uncertainty
% [      ] m = [    ]
% [Wm * H]     [Wm*h]       Wm is the weighting for the regularisation
% 
% or,    F m = f

nrays = data.ray.nrays;
nevts = data.evt.nevts;
nstas = data.stn.nstas;


%% make H
% H_damp
% [H_damp,h_damp]  = make_dampmat(data,par);
% LH_damp = svds(H_damp,1);

% H_smooth
if par.build_smooth == 0
    load H_smth H_smth;
    load H_smth LH_smth
elseif par.build_smooth == 1
   [ H_smth ] = make_smoothmat( par ); 
     LH_smth  = normest(H_smth,1e-3);
   save ('H_smth','H_smth','LH_smth','-v7.3');
end  
H_smth = [H_smth sparse(size(H_smth,1),nevts+nstas)];
LH_smth = svds(H_smth,1);

%% Data kernel

% WdGdm = (wt*ones(1,size(dGdm,2))).*dGdm;
WG = spdiags(data.ray.wt,0,nrays,nrays)*G;
Wd  = data.ray.wt.*d;
LWdGdm = svds(WG,1);

if par.scalereg
%     A = LWdGdm/LH_damp;
    B = LWdGdm/LH_smth;
else
    A = 1;
	B = 1;
end

% fprintf('norm dGdm = %.3f\n',LWdGdm);
% fprintf('norm LH_damp = %.3f\n',LH_damp);
% fprintf('norm LH_smth = %.3f\n',LH_smth);
% fprintf('scale_damp = %.3f\n',par.damp*A);
% fprintf('scale_smth = %.3f\n',par.smooth*B);
% fprintf('norm damp = %.3f\n',par.damp*A*LH_damp);
% fprintf('norm smth = %.3f\n',par.smooth*B*LH_smth);



%% Make F, f
F = [WG; par.smooth*B*H_smth];
% f = [Wres;par.damp*A*h_damp;zeros(size(H_smth,1),1)];
f = [Wd;zeros(size(H_smth,1),1)];


end

