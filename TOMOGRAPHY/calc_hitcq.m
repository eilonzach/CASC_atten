function [ par ] = calc_hitcq( data,par,K )
% [ par ] = hitmap( data,par,K )
%  function to calculate how many rays go through each voxel in the model

maxNrays = 5;
sectors  = 6;

hitc = zeros(par.nmodel,sectors);
for ii = 1:length(K.n_indx)
    sect = ceil(data.ray.baz(ii)/(360/sectors)); % find hexant (1 is 0-60, 2 is 60-120 etc.)
    hitc(K.n_indx{ii},sect) = hitc(K.n_indx{ii},sect) + 1;
end

temp = hitc;
temp(hitc>maxNrays)=maxNrays; % cap at minNrays
hitq = sum(temp,2)./(maxNrays*sectors); % quality = hits out of a possible maximum of cap*sectors

par.hitc=hitc;
par.hitq=hitq;

end

