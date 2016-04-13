function [K,data] = make_K( par, data )
fprintf('>  Making K values...  \n')

close all
tic
% stn, ray, node data
M = length(par.mx);
N = length(data.ray.d);


%------------- main loop through all rays -----------------------

sta_num = data.ray.sta_num;
slats = data.stn.lat;
slons = data.stn.lon;
selvs = data.stn.elv;
sseds = data.stn.sed;

bazs = data.ray.baz;
ray_p = data.ray.p;

Si = cell(N,1);
Sj = cell(N,1);
S  = cell(N,1);

if par.plot_raypath
rpath.rx = nan(100,N);
rpath.ry = nan(100,N);
rpath.rz = nan(100,N);
rpath.rt = nan(100,N);
end

pltry = 1;
h = waitbar(0,'Progress through ray tracing');
for iray = 1:N       % or switch to parfor to speed up
    
    % When 3-D load ray here
    
    snum = sta_num(iray); 
    % segment ends
    [rx,ry,rz,dr] = ray_trace_1D(par, ray_p(iray), bazs(iray), slats(snum), slons(snum) ,selvs(snum), sseds(snum));
    % segment midpoints
    mx = midpts(rx);
    my = midpts(ry);
    mz = midpts(rz);
    lr = diff(dr);
    
% station-individualised velocity profile
    %     rv = vel_profile(par.PS,rz,slats(snum), slons(snum) ,selvs(snum), false,sseds(snum));
    %     rt = diff(dr)./midpts(rv);
    %     sum(rt)
    mv = vel_profile(par.PS,mz+selvs(snum));
    rt = lr./mv; % travel time in each segment
    
    % function ray_nodes to find nodes potentially within RF1
    [n_indx, n_lens ] = ray_nodes_box(mx,my,mz,lr,par);
    
    if par.t_ts == 1 % if velocity inversion        
        % KERNEL VALUE IS THE PATH LENGTH 
        n_vals = n_lens;
    elseif par.t_ts == 2 % if q inversion
        % KERNEL VALUE IS THE TRAVEL TIME: LENGTH DIVIDED BY VELOCITY
        n_vals = n_lens./par.mvav(n_indx);      
    end
    
% %----------------- comment out if parallelised -----------------
%     %% PLOT RAYS
%     if par.plot_raypath
%       rpath.rx(1:length(rx),iray) = rx;
%       rpath.ry(1:length(ry),iray) = ry;
%       rpath.rz(1:length(rz),iray) = rz;
%       rpath.rt(1:length(rt),iray) = rt;
%       
%         if pltry    
%         figure(1); hold on
%         scatter3(par.mx(n_indx),par.my(n_indx),par.mz(n_indx),500*(0.000001 + n_vals),'r');
%         plot3(rx,ry,rz,'-o','LineWidth',2,'color',colour_get(data.ray.d(iray),-3,3))
%         set(gca,'zdir','reverse','xlim',[-500 500],'ylim',[-500 500],'zlim',[0 par.max_z])
%         title(sprintf('orid %u, sta %u',data.ray.orid(iray),snum))
%         
%         
%         cont = input('Press any key other than "x" to plot next ray. ','s');
%         if strcmp(cont,'x'); pltry = 0; end
%         scatter3(par.mx(n_indx),par.my(n_indx),par.mz(n_indx),500*(0.000001 + n_vals),'b');
%         end
%     end
% 
%     if mod(iray,50)==1
%     waitbar(iray/N,h)
%     end
% %----------------- comment out if parallelised -----------------
    
    Si{iray} = ones(size(n_indx))*iray;
    Sj{iray} = n_indx;
    S{iray}  = n_vals;

end


if par.plot_raypath
    L = max(sum(~isnan(rpath.rx)));
    rpath.rx(L+1:end,:)=[];
    rpath.ry(L+1:end,:)=[];
    rpath.rz(L+1:end,:)=[];
    rpath.rt(L:end,:)  =[];
    save ('rpath','rpath','-v7.3');
end

K.n_indx = Sj;
K.n_vals = S;
delete(h)

toc

