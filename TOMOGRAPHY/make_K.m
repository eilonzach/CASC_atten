function [K,data] = make_K( par, data )
fprintf('Making K values...  ')
tic
% stn, ray, node data
M = length(par.mx);
N = length(data.ray.d);


%------------- main loop through all rays -----------------------
stn_h = data.stn.h;
% ind = find(isnan(stn_h)==1);
% stn_h(ind) = 35;
sta_num = data.ray.sta_num;
gcarcs = data.ray.gcarc;
bazs = data.ray.baz;
ray_p = data.ray.p;
elats = data.ray.elat;
elons = data.ray.elon;
slats = data.stn.lat;
slons = data.stn.lon;
trdist= data.ray.trdist;
cf    = data.ray.cf;
% estimate rmax from tables
gc = gcarcs; gc(gc>90)=90; gc(gc<30)=30;
% rmax = rrmax( cf,trdist,par );


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
    if mod(iray/N,0.1)<2/N
    waitbar(iray/N,h)
    end
    
    % When 3-D load ray here
    snum = sta_num(iray); 
    [rx,ry,rz,dr] = ray_trace_1D(par.PS, ray_p(iray), bazs(iray),elats(iray),elons(iray), slats(snum), slons(snum) ,stn_h(snum), par);
%     rayxyzdr = [ rx ry rz dr ];
    rv = vel_profile(par.PS,rz,par.Vavmod,stn_h(snum));
    rmax = sqrt( (rv./cf(iray)).*dr.*(trdist(iray)-dr)./trdist(iray));
    rt = diff(dr)./midpts(rv);

    % function ray_nodes to find nodes potentially within RF1
    [n_indx, n_vals ] = ray_nodes_ze(rx,ry,rz,dr,rmax,par,cf(iray),gcarcs(iray),sum(rt));
    
    %% PLOT RAYS
    if par.plot_raypath
      rpath.rx(1:length(rx),iray) = rx;
      rpath.ry(1:length(ry),iray) = ry;
      rpath.rz(1:length(rz),iray) = rz;
      rpath.rt(1:length(rt),iray) = rt;
      
        if pltry    
        figure(1); clf; hold on
        scatter3(par.mx(n_indx),par.my(n_indx),-par.mz(n_indx),100*n_vals/mean(n_vals));
        plot3(rx,ry,-rz,'-or','LineWidth',2)
        title(sprintf('orid %u, sta %u',data.ray.orid(iray),snum))
        junk = input('Press any key other than "x" to plot next ray. ','s');
        if strcmp(junk,'x'); pltry = 0; end
        end
  
    end
    
    Si{iray} = ones(size(n_indx))*iray;
    Sj{iray} = n_indx;
    S{iray}  = n_vals;

end

% l=0;
% for k=1:length(S)
%    l=l+length(S{k});
% end
% 
% si=zeros(l,1);
% sj=zeros(l,1);
% s =zeros(l,1);
% 
% n=1;
% for k=1:length(S)
%     lng = length(Sj{k});   
%     si(n:(n+lng-1)) = Si{k};
%     sj(n:(n+lng-1)) = Sj{k};
%     s(n:(n+lng-1))  = S{k};
%     n=n+lng;
% end
% 
% G = sparse(si,sj,s,N,M);

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


