function [ model ] = model_update( model, dm, M,nevts,nstas)
% [ model ] = model_update( model, dm, M,nevts,nstas)
% 
% function to parse dm into different parts of model

model.mdq       = model.mdq + dm([1:M]);
model.estatic   = model.estatic + dm([1:nevts] + M);
model.sstatic   = model.sstatic + dm([1:nstas] + M + nevts);
end

