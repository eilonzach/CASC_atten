function [ wt ] = wtdata( data )
%[ wt ] = wtdata( data )
%   weights for data

% default: inverse of variance
wt = data.ray.sd.^(-2);
wt = wt./mean(wt);


%% FPRINTF NOTE THE IMPORTANT QC CONSTRAINTS YOU ADD HERE

% remove splits from south of -10S, where SKS tell us complex anis
% fprintf('   removing splits from south of -10S\n')
% wt(data.ray.stalat < -10 & data.ray.sect==3) = 0;

fprintf('   removing abnormally large data values\n')
wt(abs(data.ray.d) > 3.00) = 0;

wt(wt==0) = 1e-10;
end

