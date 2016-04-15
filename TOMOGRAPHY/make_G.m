function [ G ] = make_G( K,data,par )
% [ G ] = make_G( K,data,par )
%  make a matrix of derivatives of d(d)/dt, the derivatives of each datum
%  with respect to each model parm, using "model"
% 
% K is the weighting matrix, repeated twice
% 

fprintf('>  Calculating G matrix\n')

ray = data.ray;

nrays = ray.nrays;
nmod = par.nmodel;

l=0;
for k=1:length(K.n_indx)
   l=l+length(K.n_indx{k});
end


si=zeros(l,1);
sj=zeros(l,1);
s =zeros(l,1);

nn = 1;  
for iray = 1:nrays

n_indx = K.n_indx{iray};
n_vals = K.n_vals{iray};
l_nindx = length(n_indx);

inds = nn + [1:l_nindx] - 1;

si(inds,1) = iray*ones(l_nindx,1);
sj(inds,1) = n_indx;
s(inds,1)  = n_vals;
nn = nn + l_nindx;

end

del = si==0 & sj==0;
si(del) = [];
sj(del) = [];
s(del)  = [];

G = sparse(si,sj,s,nrays,nmod,l);


end

