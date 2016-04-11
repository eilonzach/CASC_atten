function [ d,G ] = add_static_terms( d,G,data,model )
% subtract event and station terms from d and G matrix
% positive e_stat or s_stat will decrease predicted tstar.
% therefore these are additional delays - positive means slow structure
% beneath just that station.
% NB adds zero rows for orids or stas are not in the data
%
% model parameter vector will be [mdq estatic sstatic]

%% event terms
rayorids = data.ray.orid;
orids = data.evt.orid;
nevts = data.evt.nevts; 

e_d = zeros(size(d));
e_G = spalloc(size(G,1),nevts,nevts);
for ie = 1:nevts
    ind = rayorids==orids(ie);
    e_d(ind) = model.estatic(ie);
    e_G(:,ie) = -ind;
end

%% station terms
raystas = data.ray.sta_num;
stas = data.stn.num;
nstas = data.stn.nstas;

s_d = zeros(size(d));
s_G = spalloc(size(G,1),nstas,nstas);
for is = 1:nstas
    ind = raystas==stas(is);
    s_d(ind) = model.sstatic(is);
    s_G(:,is) = -ind;
end

d = d - e_d - s_d; % positive e_stat or s_stat will decrease tstar_pred
G = [G,e_G,s_G];


end

