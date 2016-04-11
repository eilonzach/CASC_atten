function [ plot_model ] = conv2plotable( model,par )
%[ plot_model ] = conv2plotable( model,par )
%   convert the model in vx,vy,vz to vv,va,vp and parse into shape for
%   plotting
shape = [par.ny,par.nx,par.nz];

plt.x = reshape(par.mx,shape);
plt.y = reshape(par.my,shape);
plt.z = reshape(par.mz,shape);
plt.lt= reshape(par.mlt,shape);
plt.ln= reshape(par.mln,shape);

if par.t_ts == 1
    plt.val = reshape(model.mdv,shape);
elseif par.t_ts ==2;
    plt.val = reshape(model.mdq,shape);
end

%% hit quality
try
    plt.hq = reshape(par.hitq,shape);
catch
    plt.hq = nan(size(plt.val));
end

%% semblance
try 
    plt.semb = reshape(par.semb,shape);
catch
    plt.semb = nan(size(plt.val));
end

plot_model = plt;
end

