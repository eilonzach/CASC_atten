clear all
close all
global F f G res
% profile on
cd ~/Documents/MATLAB/CASC_atten/TOMOGRAPHY


fprintf('\n=========== RUNNING INVERSION ===========\n\n')

%% ALL PARMS ESTABLISHED IN PARMFILE
% consider moving frequently changed ones to this script
    fprintf('>  Establishing parameters\n')
run('PARMS')

%% read data
    fprintf('>  Reading in data and station details\n')
data = read_data(datfile,stafile,par);

%% weight data
if par.wtdata
    fprintf('Weighting of data being altered by function wtdata\n')
    fprintf('CHECK wtdata CAREFULLY TO SEE WHAT IT IS DOING!!\n')    
    data.ray.wt = wtdata(data,par);
else
    fprintf('>  No data weighting\n')
    data.ray.wt = ones(size(data.ray.d));
end

%% delete bad orids!
% [ data ] = wipeorids( data );

%% setup tomo geom
    fprintf('>  Setting up geometry\n')
[ data,par ] = setup_geom( data,par );

%% Crustal correction
% if par.crust_corr
%     data = crust_corr(data,par);
% end

%% Plot data
if par.plot_data
    fprintf('>  Plotting data\n')
    plot_data(data,par,1)
end


%% Make starting model - required to put starting model into parm
[model,par] = make_start_model(par,data);

model_1 = model;

if par.plot_inmodel && ~par.synth_test
    [plot_inmodel] = conv2plotable(model,par);
    if par.force2D==0
        plot_basic(plot_inmodel,par,2,par.saveopt)
    else
        plot_zslice(plot_inmodel,par,2,par.saveopt)
    end
end

%% =========== Make Kernel =============================================
if par.t_ts == 1
    Kfile = 'K_v';
elseif par.t_ts == 2
    Kfile = 'K_q';
end

if par.build_K == 0
    load(Kfile);
    if length(K.n_indx)~=data.ray.nrays || max(K.n_indx{1}) > par.nmodel
        error('Loaded K not appropriate for this data... build K again')
    end
elseif par.build_K == 1
    [K,data] = make_K( par, data ); 
    save (Kfile,'K','-v7.3');
end



%% Hit quality?
par = calc_hitcq( data,par,K );
if par.plot_hitq
    plot_hitq(par);
end

%% Raypaths
if par.plot_raypath
    plot_raypaths(par,data,1);
    plot_raypaths(par,data,1);
end

%% Synthetic test
if par.synth_test
    synth_test
    return
end


%% ==================== DO INVERSION ================= %%

G = make_G( K,data,par );

%% predicted data from starting model
d_pred = G*model.mval;
[d_pred,G] = add_static_terms(d_pred,G,data,model);

%% calculate residual and do inverse problem
d_use = data.ray.d;
res = d_use - d_pred;

%% ADD REGULARISATION
[ F,f ] = make_F_f( G,res,data,par );

%% SOLVE
if strcmp(par.solver,'lsqr')
    fprintf('>  Solving F*dm = f using LSQR\n')
    [dm,flag,relres,iter,resvec] = lsqr( F, f, 1e-5, 500 );
    if iter==500, fprintf('Warning LSQR not converging\n'); end
elseif strcmp(par.solver,'bicg')
    fprintf('>  Solving F*dm = f using biconjugate gradient method\n')
    [ dm ] = solve_bicg( F, f, 1e-4, 1e4 );
end

%% UPDATE
[ model ] = model_update( model, dm, par.nmodel,data.evt.nevts,data.stn.nstas);

%% Final residual
d_pred = G*[model.mval;model.estatic;model.sstatic];
res = d_use - d_pred;

fprintf('>  Results summary:\n')
% profile viewer


%% RESULT STATS
vr  = variance_reduction(d_use,d_pred);
wvr = variance_reduction(data.ray.wt.*d_use,data.ray.wt.*d_pred);

fprintf('Final RMS error: %.2f\n',rms(res));
fprintf('Variance reduction = %.2f %%\n',vr);
fprintf('Weighted variance reduction = %.2f %%\n',wvr);

fprintf('RMS event static = %.2f s \n',rms(model.estatic))
fprintf('RMS station static = %.2f s \n',rms(model.sstatic))

plot_results_info(par,model,data,d_use,res)
plot_model = conv2plotable(model,par);
plot_basic(plot_model,par,1,par.saveopt)
return
% [model_hq ] = hiQmodel( model,model_1,par,0.3 );
% d_pred_hq = G*[model_hq.mval;model_hq.estatic;model_hq.sstatic];
% res_hq     = d_use - d_pred_hq;
% vr_hq = variance_reduction(d_use,d_pred_hq);
% 
% fprintf('HI-Q Variance reduction = %.2f %%\n',vr_hq);

% some squeezing test stats...
fprintf('\nSqueezing analysis:\n')
if par.niter>0
maxitersq = ceil(par.squeeze*par.niter);
fprintf('With squeezing to %.0f km, achieve\n', par.zsqz)
fprintf('  %.2f %% of final total var reduction\n',100*vr.all(maxitersq)/vr.all(end))
fprintf('  %.2f %% of final tdiff var reduction\n',100*vr.tdiff(maxitersq)/vr.tdiff(end))
fprintf('  %.2f %% of final dT var reduction\n',100*vr.dT(maxitersq)/vr.dT(end))
fprintf('  %.2f %% of final wtvar reduction\n',100*wvr.all(maxitersq)/wvr.all(end))
end

indsqz = par.mz>par.zsqz; 
fprintf('Norm of dV structure below zsqz = %.2f ',norm(model.dv(indsqz)));
fprintf(' = %.2f%% of a total norm of %.2f\n',100*(norm(model.dv(indsqz))./norm(model.dv))^2,norm(model.dv));
fprintf('Norm of A structure below zsqz = %.2f ',norm(model.aa(indsqz)));
fprintf(' = %.2f%% of a total norm of %.2f\n',100*(norm(model.aa(indsqz))./norm(model.aa))^2,norm(model.aa));

%% plots etc.
load('semblance');
par.semb = semb;

save ('model','model','-v7.3');
save data data; 
save par par
if par.force2D==1
    save('model_2D','model','-v7.3')
    save par_2D par
end

if par.plot_outputs

plot_results_info(par,model,dm,data,d_use,res)
plot_model = conv2plotable(model,par);
if par.force2D==1
    plot_Zmaps(plot_model,par,data,par.saveopt)
else
    plot_Hmaps(plot_model,par,1,par.saveopt)
    if par.plot_Zslices
    plot_Zmaps(plot_model,par,data,par.saveopt)
    end
end

end


%% Manual optional extra plots
return
plot_errors( data,par,res )
plot_pretty(plot_model,par,1,1)
plot_zslice(plot_model,par,1,1)
plot_xy(plot_model,par,1,1)

return

