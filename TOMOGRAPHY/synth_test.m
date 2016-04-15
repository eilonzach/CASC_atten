fprintf('\nSTARTING SYNTHETIC TEST\n')
% profile clear
% profile on
%% ================  setup synthetic model  ================  
synth_model = make_synth_model(par,data);

if par.plot_inmodel
    [plot_simodel] = conv2plotable(synth_model,par);
    plot_basic(plot_simodel,par,2,par.saveopt)
end

mf = synth_model.mval;

%% Calc synthetic data
fprintf('>  Calculating synthetic data\n')
% data kernel
G = make_G( K,data,par );
% make synthetic data
d_synth = G*mf;
% event and station terms
% noise
if par.synth_noisy
    d_synth = d_synth + random('norm',0,par.sym.noise,size(d_synth));
end
[d_synth,~] = add_static_terms(d_synth,G,data,synth_model);

%% Plot data
if par.plot_data
    d_real = data.ray.d; % save the real data
    data.ray.d = d_synth; % put synthetic data into ray struct to plot
    fprintf('>  Plotting data\n')
    plot_data(data,par,1) % plot synthetic data
    data.ray.d = d_real; % put real data back in
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
[d_pred,G] = add_static_terms(d_pred,G,data,model);

%% calculate residual and do inverse problem
d_use = d_synth;
res = d_use - d_pred;

%% ADD REGULARISATION
[ F,f ] = make_F_f( G,res,data,par );

%% SOLVE
if strcmp(par.solver,'lsqr')
    fprintf('>  Solving F*dm = f using LSQR\n')
    [dm,flag,relres,iter,resvec] = lsqr( F, f, 1e-5, 1000 );
    if iter==1000, fprintf('Warning LSQR not converging\n'); end
elseif strcmp(par.solver,'bicg')
    fprintf('>  Solving F*dm = f using biconjugate gradient method\n')
    [ dm ] = solve_bicg( F, f, 1e-4, 1e4 );
end

%% UPDATE
[ model ] = model_update( model, dm, par.nmodel,data.evt.nevts,data.stn.nstas);

toc
solved = true;
end

par.si_model = synth_model;
par.so_model = model;
% par = zelt_semblance_anis(par);

%% Final residual
d_pred = G*[model.mval;model.estatic;model.sstatic];
res = d_use - d_pred;


%% ================== Plotting ==================

plot_results_info(par,model,data,d_use,res,synth_model)
% profile viewer
return
if par.plot_synout
    [plot_somodel] = conv2plotable(model,par);
    plot_basic(plot_somodel,par,3,par.saveopt)
    plot_Zmaps(plot_model,par)
end
return

if strcmp(par.sym.opt,'checker')
    semb = par.semb;
    save('semblance','semb')
end

save par_synth par

return
