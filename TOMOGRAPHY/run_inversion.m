% clear all
% close all
global F f dGdm res
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
% model = ideal_model;

if par.plot_inmodel && ~par.synth_test
    [plot_inmodel] = conv2plotable(model,par);
    if par.force2D==0
        plot_basic(plot_inmodel,par,2,par.saveopt)
    else
        plot_zslice(plot_inmodel,par,2,par.saveopt)
    end
end

% % special starting values
% model.sstatic(strcmp(data.stn.sta,'J')) = 1; % slow sediments?
% model.sstatic(strcmp(data.stn.sta,'H')) = 1; % slow sediments?

%% =========== Make Kernel =============================================
if par.build_K == 0
    load K K;
    if length(K.n_indx)~=data.ray.nrays
        error('Loaded K not appropriate for this data... build K again')
    end
elseif par.build_K == 1
    [K,data] = make_K( par, data ); 
    save ('K','K','-v7.3');
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
return


%% ==================== DO INVERSION ================= %%
allres = zeros(par.niter,1);
wtres  = zeros(par.niter,1);
tabsres= zeros(par.niter,1);
dtsres = zeros(par.niter,1);
vr = struct('all',zeros(par.niter,1),'tdiff',zeros(par.niter,1),'dT',zeros(par.niter,1));
wvr = struct('all',zeros(par.niter,1),'tdiff',zeros(par.niter,1),'dT',zeros(par.niter,1));

for kk = 1:par.niter
par.iter=kk;
fprintf('\nIteration number %u\n',kk)
tic
%% Predicted data
gm   = make_gm( K,data,model,par ); 
dGdm = make_dGdm( K,data,model,par );

if any(imag(gm)~=0) | any(imag(dGdm)~=0) %#ok<OR2>
    fprintf('Imaginary data!!')
    return
end
    
%% event and station terms
[gm,dGdm] = add_static_terms(gm,dGdm,data,model);

%% calculate residual and do inverse problem
if par.crust_corr
    d_use = data.ray.d - data.ray.corr1;
else
    d_use = data.ray.d;
end
res = d_use - gm;
allres(kk) = norm(res);
wtres(kk) = norm(data.ray.wt.*res);
tabsres(kk) = norm(res(data.ray.sect<3));
dtsres(kk) = norm(res(data.ray.sect==3));
if kk>1
    vr.all(kk-1) = variance_reduction(d_use,gm);
    vr.tdiff(kk-1) = variance_reduction(d_use(data.ray.sect<3),gm(data.ray.sect<3));
    vr.dT(kk-1) = variance_reduction(d_use(data.ray.sect==3),gm(data.ray.sect==3));
    wvr.all(kk-1) = variance_reduction(data.ray.wt.*d_use,data.ray.wt.*gm);
    wvr.tdiff(kk-1) = variance_reduction(data.ray.wt(data.ray.sect<3).*d_use(data.ray.sect<3),...
                                         data.ray.wt(data.ray.sect<3).*gm(data.ray.sect<3));
    wvr.dT(kk-1) = variance_reduction(data.ray.wt(data.ray.sect==3).*d_use(data.ray.sect==3),...
                                         data.ray.wt(data.ray.sect==3).*gm(data.ray.sect==3));
    % fprintf('Variance reduction = %.2f %%\n',100*(allres(1).^2 - allres(kk).^2)/allres(1)^2);
    % fprintf('Weighted variance reduction = %.2f %%\n',100*(wtres(1).^2 - wtres(kk).^2)/wtres(1)^2);
    fprintf('After iter %.0f Variance reduction = %.2f %%\n',kk-1,vr.all(kk-1));
    fprintf('After iter %.0f Weighted variance reduction = %.2f %%\n',kk-1,wvr.all(kk-1));
end

% END if norm res is small enough
figure(1);
plot(allres,'o-k');
hold on
plot(tabsres,'o-r');
plot(dtsres,'o-b');
hold off
title('k=all,r=tabs,b=dts')
xlim([1 max([kk+1,10])])
pause(0.001)

if kk>1
    if abs(allres(kk)-allres(kk-1))/abs(allres(kk-1)) < par.stopiter
        fprintf('\n NO CHANGE IN RES... STOPPING\n')
        break
    end
end

%% ADD REGULARISATION
[ F,f ] = make_F_f( dGdm,res,data,par );

%% SOLVE
if strcmp(par.solver,'lsqr')
    fprintf('Solving F*dm = f using LSQR\n')
    [dm,flag,relres,iter,resvec] = lsqr( F, f, 1e-6, 500 );
    if iter==500, fprintf('Warning LSQR may not be converging\n'); end
elseif strcmp(par.solver,'bicg')
    fprintf('Solving F*dm = f using biconjugate gradient method\n')
    [ dm ] = solve_bicg( F, f, 1e-4, 1e4 );
end


%% UPDATE
[ model ] = model_update( model, dm, par.nmodel,data.evt.nevts,data.stn.nstas);

%% Plot changes
if par.plot_everyNtime > 0
if mod(kk-1,par.plot_everyNtime)==0
    plot_results_info(par,model,dm,data,d_use,res)
    
    [plot_model] = conv2plotable(model,par);
    plot_basic(plot_model,par,1,0)
    
    pause(0.01)
end
end

toc
par.build_K         = 0;
par.build_smooth    = 0;

end % loop on inversion iterations
fprintf('\nResults summary:\n')
% profile viewer

[model.vv,model.aa] = va2xy(model.mvx,model.mvy,'backward');
model.dv = 100*(model.vv-par.mvav)./par.mvav;

%% RESULT STATS
gm      = make_gm( K,data,model,par ); 
[gm,~] = add_static_terms(gm,dGdm,data,model);
res     = d_use - gm;
vr.all(end)   = variance_reduction(d_use,gm);
vr.tdiff(end)   = variance_reduction(d_use(data.ray.sect<3),gm(data.ray.sect<3));
vr.dT(end)   = variance_reduction(d_use(data.ray.sect==3),gm(data.ray.sect==3));
wvr.all(end) = variance_reduction(data.ray.wt.*d_use,data.ray.wt.*gm);
wvr.tdiff(end) = variance_reduction(data.ray.wt(data.ray.sect<3).*d_use(data.ray.sect<3),...
                                     data.ray.wt(data.ray.sect<3).*gm(data.ray.sect<3));
wvr.dT(end) = variance_reduction(data.ray.wt(data.ray.sect==3).*d_use(data.ray.sect==3),...
                                         data.ray.wt(data.ray.sect==3).*gm(data.ray.sect==3));
fprintf('Final RMS error: %.2f\n',rms(res));
% fprintf('Variance reduction = %.2f %%\n',100*(allres(1).^2 - norm(res).^2)/allres(1)^2)
% fprintf('Weighted variance reduction = %.2f %%\n',100*(wtres(1).^2 - norm(data.ray.wt.*res).^2)/wtres(1)^2)
% fprintf('Only tdiff variance reduction = %.2f %%\n',100*(tabsres(1).^2 - norm(res(data.ray.sect<3)).^2)/tabsres(1)^2)
% fprintf('Only dt variance reduction = %.2f %%\n',100*(dtsres(1).^2 - norm(res(data.ray.sect==3)).^2)/dtsres(1)^2)

fprintf('Variance reduction = %.2f %%\n',vr.all(end));
fprintf('Weighted variance reduction = %.2f %%\n',wvr.all(end));
fprintf('Only tdiff variance reduction = %.2f %%\n',vr.tdiff(end));
fprintf('Only dt variance reduction = %.2f %%\n',vr.dT(end));

fprintf('RMS station diff-time static = %.2f s \n',0.5*rms(model.sstatic(1:data.stn.nstas) + model.sstatic(data.stn.nstas+1:end)))
fprintf('RMS station splitting static = %.2f s \n',rms(model.sstatic(1:data.stn.nstas) - model.sstatic(data.stn.nstas+1:end)))


[model_hq ] = hiQmodel( model,model_1,par,0.3 );
gm_hq   = make_gm( K,data,model_hq,par ); 
[gm_hq,~] = add_static_terms(gm_hq,dGdm,data,model);
res_hq     = d_use - gm_hq;
vr_hq = variance_reduction(d_use,gm_hq);

% fprintf('HI-Q Variance reduction = %.2f %%\n',100*(allres(1).^2 - norm(res_hq).^2)/allres(1)^2)
fprintf('HI-Q Variance reduction = %.2f %%\n',vr_hq);

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
    plot_zslice(plot_model,par,1,par.saveopt)
else
    plot_basic(plot_model,par,1,par.saveopt)
    if par.plot_Zslices
    plot_zslice(plot_model,par,1,par.saveopt)
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

