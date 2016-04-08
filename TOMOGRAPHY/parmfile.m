% clear all
close all
cd('/Users/zeilon/Documents/MATLAB/CASC_atten/TOMOGRAPHY');
%% Set parms
fprintf('Running parameters file')
%% data files
datfile = 'data/data_dtstar_ST.dat';
stafile  = 'data/stations.dat';
 
par=struct([]);
par(1).PS = 2;         % 1 for P, 2 for S
 
%% plotting parms - zero for speed
par.plot_crustcorr  = 0;
par.plot_data       = 1;
par.plot_inmodel    = 0;
par.plot_synout     = 1;
par.plot_Zslices    = 0;
par.plot_hitq       = 0;
par.plot_raypath    = 0;
par.plot_outputs    = 1;
par.plot_everyNtime = 0; % N>0 will plot each N iterations
 
%% inversion running parms - zero for speed
REBUILD             = 0;
par.redo_crust_corr = REBUILD; % <== REDO IF DATA/MODEL CHANGES
par.build_K         = REBUILD; % <== REDO IF DATA/MODEL CHANGES
par.build_smooth    = 1; % <== REDO IF DATA/MODEL CHANGES
par.solver          = 'lsqr'; % 'lsqr' or 'bicg'

par.synth_test      = 0;  % this will stop actual inversion
par.synth_noisy     = 1;  % option to include 'realistic' noise in synthetic data

par.niter           = 10; %5
par.stopiter        = 0.00; %if fractional res reduction does not exceed this, cease
par.crust_corr      = 0; % option to do crustal correction
par.wtdata          = 1; % extra QCs found in wtdata.m function

par.force2D         = 0.1; %0.2% 0=no, 1=yes, 0<frac<1=force 2D for first frac of runs

par.squeeze         = 0.2; %0.2% 0=no, 1=yes, 0<frac<1=force 2D for first frac of runs
par.zsqz            = 175; %180% squeezing depth in km, if -ive then sqz to below this depth

par.saveopt         = 1;

% regularisation parms
par.damp            = 3; %3% scaling for whole damp mat: more damping parms found in dampmat.m
par.smooth          = 3;  %3%  scaling for whole smth mat: good result with 1 - 

par.scalereg        = 0; %       option to scale regularisation matrices in make_F_f
% ==== more damping parms in make_dampmat.m ===%
 
%% model space parms
    par.origin  = [44 -120];   % model origin
    mstruct = defaultm('mercator');
    mstruct.origin = [par.origin 0];
    mstruct = defaultm( mstruct ); 
    par.map_proj = mstruct;
 
    par.dh      = 30; %30   % horizontal spaces for nodes
    par.dh_max  = 80; %80   % horizontal spacing to which it scales, at twice array aperture

    par.dz      = 20; %35   % vertical spacing for nodes beneath moho
    par.dz_max  = 35; %40   % vertical spacing to which it scales, at model base
 
    par.zmax    = 400; %300 % maximum vertical distance
    par.zmin    = 30;  %40  % top of model (above this must be accounted for by crust
 
% kernel integration calculation parms
    par.kidd    = 10;      % along-ray segment lengths to calc at
    par.kidr    = 1;       % ray-normal distances to calc at
 
%% starting model parms
    par.Asstart = 0; % percentage anisotropy for starting model
    par.Vavmod  = 2; % starting model: 1=AK135,2=Raj+Ferris,3=STW105
    
%% synth model parms
    % custom structure
    sym.acx = [  0	 70     -50]           ;% vector of x-coords for anomaly centres
    sym.acy = [-40	  0     100]           ;% vector of y-coords for anomaly centres
    sym.acz = [ 110	 90      90]           ;% vector of z-coords for anomaly centres
    
    sym.awx = [250	 60     180]           ;% vector of x-widths for anomalies (km)
    sym.awy = [ 80	180      80]           ;% vector of y-widths for anomalies (km)
    sym.awz = [100	130      70]           ;% vector of z-widths for anomalies (km)

    sym.adv = [ -6	  0       4]           ;% vector of v pertubations of anomalies (percent)
    sym.aa  = [  5	  0      -4]           ;% vector of anisotropy values of anomalies (percent)
    
    % checker structure
    sym.dv = 5;
    sym.da = 4;
    sym.stag = 1; % 1 or 0 - option to stagger anis and vel checkers (if checker model)
    
    sym.opt  = 'checker'; % 'checker' or 'custom'
    
    sym.noise = 0.1;                     % 0.05 % std. of gaussian perturbations to data values (s)
    
    par.sym = sym;
    
    
    
    
%     % if 2D and checker, make 2D checkers...
%     if par.force2D && strcmp(sym.opt,'checker')
% 	sym.acx = [   0    0    0    0    0    0    0    0     0    0    0    0    0    0 ]           ;% vector of x-coords for anomaly centres
%     sym.acy = [-235 -235  -60  -60   60   60  200  200  -130 -130    0    0  120  120 ]           ;% vector of y-coords for anomaly centres
%     sym.acz = [  75  155   75  155   75  155   75  155    75  155   75  155   75  155 ]           ;% vector of z-coords for anomaly centres
%     
%     sym.awx = [ 1e3  1e3  1e3  1e3  1e3  1e3  1e3  1e3   1e3  1e3  1e3  1e3  1e3  1e3 ]           ;% vector of x-widths for anomalies (km)
%     sym.awy = [ 150  150   90   90   90   90  120  120    90   90   90   90   90   90 ]           ;% vector of y-widths for anomalies (km)
%     sym.awz = [  50   50   50   50   50   50   50   50    50   50   50   50   50   50 ]           ;% vector of z-widths for anomalies (km)
% 
%     sym.adv = [  -5    5    5   -5   -5    5    5   -5     0    0    0    0    0    0 ]           ;% vector of v pertubations of anomalies (percent)
%     sym.aa  = [   0    0    0    0    0    0    0    0    -5    5    5   -5   -5    5 ]           ;% vector of anisotropy values of anomalies (percent)
%     sym.opt  = 'custom'; % 'checker' or 'custom'
%     
%     par.sym = sym;
%     end
    
% %     
% %% RUN!
% if ~exist('bokachoda','var')
% run_inversion
% end
% 
% % 
