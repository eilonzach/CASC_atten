fprintf('\nSTARTING SYNTHETIC TEST\n')

%% setup synthetic model
synth_model = make_synth_model(par);

if par.plot_inmodel
    [plot_simodel] = conv2plotable(synth_model,par);
    if par.force2D==1
        plot_zslice(plot_simodel,par,2,par.saveopt)
    else
        plot_basic(plot_simodel,par,2,par.saveopt)
    end
end

if par.t_ts == 1; 
    mf = synth_model.mdv;
elseif par.t_ts == 2; 
    mf = synth_model.mdq;
end

%% data kernel
G = make_G( K,data,par );

d_obs = G*mf;

[d_obs,G] = add_static_terms(d_obs,G,data,model_1);


if par.synth_noisy
    d_obs = d_obs + random('norm',0,par.sym.noise,size(d_obs));
end

%% ==================== DO INVERSION ================= %%
%% ==================== DO INVERSION ================= %%
[model_1,par] = make_start_model(par,data);
% Tweak parms for synth case
% par.damp_evt = 1e4;
% par.damp_stn = 1e4;
% par.damp_evt_anis = 1e4;
% par.damp_stn_anis = 1e4;
mf = model_1.mdq;

allres = zeros(par.niter,1);

% ------------------------- START ITERATIONS ---------------------%

for kk = 1:par.niter
par.iter=kk;
fprintf('\nIteration number %u\n',kk)
tic

%% Predicted data
d = G*[model_1.mdq;model_1.estatic;model_1.sstatic];

% if any(imag(d)~=0) | any(imag(G)~=0) %#ok<OR2>
%     fprintf('Imaginary data!!')
%     return
% end

%% event and station terms

%% calculate residual and do inverse problem
d_use = d_obs;
res = d_use - d;
allres(kk) = norm(res);

fprintf('Variance reduction = %.2f %%\n',variance_reduction(d_use,d));

% END if norm res is small enough
figure(1);clf;plot(allres,'o-'); pause(0.001)
if kk>1
    if abs(allres(kk)-allres(kk-1)) < 0.01
        fprintf('\n NO CHANGE IN RES... STOPPING\n')
        break
    end
end

%% ADD REGULARISATION
[ F,f ] = make_F_f( G,res,data,par );

%% SOLVE
fprintf('Solving F*dm = f using LSQR\n')
[dm,flag,relres,iter,resvec] = lsqr( F, f, 1e-4, 400 );
if iter==400, fprintf('Warning LSQR may not be converging\n'); end
%% UPDATE
[ model_1 ] = model_update( model_1, dm, par.nmodel,data.evt.nevts,data.stn.nstas);

%% Plot changes
if par.plot_everyNtime > 0
if mod(kk-1,par.plot_everyNtime)==0
    plot_results_info(par,model_1,dm,data,d_use,res)
    
    [plot_somodel] = conv2plotable(model_1,par);
    plot_basic(plot_somodel,par,3,0)
    
    pause(0.01)
end
end

toc
par.build_K         = 0;
par.build_smooth    = 0;

end % loop on inversion iterations

par.si_model = synth_model;
par.so_model = model_1;
par = zelt_semblance_anis(par);

plot_results_info(par,model_1,dm,data,d_use,res)

if par.plot_synout
    [plot_somodel] = conv2plotable(model_1,par);
    if par.force2D==1
        plot_zslice(plot_somodel,par,3,par.saveopt)
    else
        plot_basic(plot_somodel,par,3,par.saveopt)
    end
end

if strcmp(par.sym.opt,'checker')
    semb = par.semb;
    save('semblance','semb')
end

save par_synth par

return
