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


%% H_damp
[H_damp,h_damp]  = make_dampmat(data,par);

%% H_smooth
if par.t_ts == 1
    Hsmthfile = 'H_smth_v';
elseif par.t_ts == 2
    Hsmthfile = 'H_smth_q';
end

if par.build_smooth == 0
    load(Hsmthfile);
elseif par.build_smooth == 1
   [ H_smth ] = make_smoothmat( par ); 
   save (Hsmthfile,'H_smth','-v7.3');
end  
H_smth = [H_smth sparse(size(H_smth,1),nevts+nstas)]; % add station+evt columns

%% Data kernel
WG = spdiags(data.ray.wt,0,nrays,nrays)*G;
Wd  = data.ray.wt.*d;

% scale regularisation
if par.scalereg
    LWdGdm = svds(WG,1);
    LH_damp = svds(H_damp,1);
    LH_smth  = normest(H_smth,1e-3);
% LH_smth = svds(H_smth,1);
    
    A = LWdGdm/LH_damp;
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
F = [WG; A*H_damp; B*par.smooth*H_smth];
f = [Wd; A*h_damp; zeros(size(H_smth,1),1)];

% pt_k = [45,-125,80]
% [ r_k ] = Rmatrix( F, par, pt_k )
% pause
end

