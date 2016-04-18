clear all
close all
global F f G res
% profile on
cd ~/Documents/MATLAB/CASC_atten/TOMOGRAPHY

%% L test parms
res_redo = 0;
saveopt = 0;
%penalty function
resid2rough = 10; %relative scaling in penalty between roughness and residual 
% (higher minimises residual, lower minimises roughness)
% preferred value around 0.5 (since resids are higher already)
odir = '~/Documents/MATLAB/CASC_atten/TOMOGRAPHY/results/';

Ltest = struct([]);

Ltest(1).damp = 3*logspace(-2,1,9)'; % was logspace(-2,2,21)
Ltest.smooth =  [.1,1,3,10,20,30]; % was logspace(-2,2,21)

Nd = length(Ltest.damp);
Ns = length(Ltest.smooth);

Ltest.resid     = nan(Nd,Ns);
Ltest.norm      = nan(Nd,Ns);
Ltest.vr        = nan(Nd,Ns);
Ltest.vr_hq     = nan(Nd,Ns);
Ltest.wvr       = nan(Nd,Ns);
Ltest.norm_estt = nan(Nd,Ns);
Ltest.norm_sstt = nan(Nd,Ns);
Ltest.norm_dat  = nan(Nd,Ns);

fprintf('\n=========== RUNNING L-CURVE INVERSION ===========\n\n')

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
    plot_Hslice(plot_inmodel,par,2,par.saveopt)
    plot_zslice(plot_inmodel,par,2,par.saveopt)
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
end

%% ==================== DO INVERSION ================= %%

G = make_G( K,data,par );

%% predicted data from starting model
d_pred = G*model.mval;
[d_pred,G] = add_static_terms(d_pred,G,data,model);

%% calculate residual and do inverse problem
d_use = data.ray.d;
res = d_use - d_pred;


%% ---------------------  LTESTING  ---------------------
%% ---------------------  LTESTING  ---------------------
fprintf('\n>  Starting L-curve iterations\n')
Ltest.par = par;
for id = 1:length(Ltest.damp)
for is = 1:length(Ltest.smooth)
    par.damp = Ltest.damp(id);
    par.smooth = Ltest.smooth(is);
    model = model_1;
    d_pred = G*[model.mval;model.estatic;model.sstatic];
    res = d_use - d_pred;
    fprintf('    Damp=%.2f  Smooth=%.2f\n',par.damp,par.smooth)
    
    %% ADD REGULARISATION
    [ F,f ] = make_F_f( G,res,data,par );

    %% SOLVE
    if strcmp(par.solver,'lsqr')
        [dm,flag,relres,iter,resvec] = lsqr( F, f, 1e-5, 500 );
%         if iter==500, fprintf('Warning LSQR not converging\n'); end
    elseif strcmp(par.solver,'bicg')
%         fprintf('>  Solving F*dm = f using biconjugate gradient method\n')
        [ dm ] = solve_bicg( F, f, 1e-4, 1e4 );
    end

    %% UPDATE
    [ model ] = model_update( model, dm, par.nmodel,data.evt.nevts,data.stn.nstas);

    %% Final residual
    d_pred = G*[model.mval;model.estatic;model.sstatic];
    res = d_use - d_pred;

    %% RESULT STATS
    vr  = variance_reduction(d_use,d_pred);
    wvr = variance_reduction(data.ray.wt.*d_use,data.ray.wt.*d_pred);

    %% Put into Ltest struct
    Ltest.resid(id,is)     = norm(res);
    Ltest.norm(id,is)      = norm(model.mval);
    Ltest.vr(id,is)        = vr;
    Ltest.vr_hq(id,is)     = nan;
    Ltest.wvr(id,is)       = wvr;
    Ltest.norm_estt(id,is) = norm(model.estatic);
    Ltest.norm_sstt(id,is) = norm(model.sstatic);
    Ltest.norm_dat(id,is)  = norm(d_pred);

    fprintf('      Variance reduction: %.2f\n',Ltest.vr(id,is))
    fprintf('      Misfit norm:        %.2f\n',Ltest.resid(id,is))
    fprintf('      Model norm:         %.2f\n',Ltest.norm(id,is))
end  
end

datstr = {'dT','dtstar'};
PSstr = {'PZ','ST'};
ofile = [odir,'Ltest_',datstr{par.t_ts},'_',PSstr{par.PS}];

save(ofile,'Ltest')

plot_Lcurves

