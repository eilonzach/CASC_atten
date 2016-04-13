fprintf('\nSTARTING SYNTHETIC TEST\n')

%% ================  setup synthetic model  ================  
synth_model = make_synth_model(par);

if par.plot_inmodel
    [plot_simodel] = conv2plotable(synth_model,par);
    plot_basic(plot_simodel,par,2,par.saveopt)
end

mf = synth_model.mval;

%% Calc synthetic data
% data kernel
G = make_G( K,data,par );
% make synthetic data
d_synth = G*mf;
% event and station terms
% noise
if par.synth_noisy
    d_synth = d_synth + random('norm',0,par.sym.noise,size(d_synth));
end

%% ==================== DO INVERSION ================== %%
%% ============ copy from here for real thing ========= %%

[model,par] = make_start_model(par,data);

%% data kernel (redundant, but want to have the full inversion structure in here
G = make_G( K,data,par );

tic
count = 0;
solved = false;
while solved == false % loop if you will do squeezing
count = count+1;

%% Predicted data
d_pred = G*model.mval;
[d_synth,G] = add_static_terms(d_synth,G,data,model);

%% calculate residual and do inverse problem
d_use = d_synth;
res = d_use - d_pred;

%% ADD REGULARISATION
[ F,f ] = make_F_f( G,res,data,par );

%% SOLVE
fprintf('Solving F*dm = f using LSQR\n')
[dm,flag,relres,iter,resvec] = lsqr( F, f, 1e-4, 400 );
if iter==400, fprintf('Warning LSQR may not be converging\n'); end

%% UPDATE
[ model ] = model_update( model, dm, par.nmodel,data.evt.nevts,data.stn.nstas);

toc
solved = true;
end

par.si_model = synth_model;
par.so_model = model;
% par = zelt_semblance_anis(par);

%% ================== Plotting ==================

% plot_results_info(par,model_1,dm,data,d_use,res)


if par.plot_synout
    [plot_somodel] = conv2plotable(model,par);
    plot_basic(plot_somodel,par,3,par.saveopt)
end
return

if strcmp(par.sym.opt,'checker')
    semb = par.semb;
    save('semblance','semb')
end

save par_synth par

return
