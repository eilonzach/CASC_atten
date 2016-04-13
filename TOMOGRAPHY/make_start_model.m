function [ model, par ] = make_start_model( par,data)
%[ model, par ] = make_start_model( par,nevts,nstas  )
%  make the starting model for the inversion


nevts = data.evt.nevts;
nstas = data.stn.nstas;

model = struct('mval',zeros(par.nmodel,1));
           
% positive e_stat or s_stat will decrease predicted time.
% therefore these are additional delays - positive means slow structure
% beneath just that station.
           
% static terms
model.estatic = zeros(nevts,1); 
model.sstatic = zeros(nstas,1); 

par.start_model = model;
end

