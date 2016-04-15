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

if par.t_ts == 1;
    % calculated perturbations in slowness, but want perturabations in V
    model.mval = dQ_to_dq(model.mval); % same conversion scheme as from dq to dQ
    plt.val = reshape(model.mval,shape);
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

