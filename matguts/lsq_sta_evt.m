function [ sta_terms,evt_terms ] = lsq_sta_evt( all_mat, epsilon )
% Script to solve the damped least squares problem for best fitting station and
% event values from some matrix of values that are (at each station) a
% combination of station and event terms.
% 
%  Usage:
%       [ sta_terms,evt_terms ] = lsq_sta_evt( all_mat )
% 
%  Inputs:
%       all_mat   = nstas x nevts matrix of values (possibly with nans)
%       epsilon   = strength of damping for event terms.
% 
%  Outputs:
%       sta_terms = nstas x 1 vector of best-fitting station values
%       evt_terms = nevts x 1 vector of best-fitting event values

if nargin < 2
    epsilon = 0.001;
end

[nstas,nevts] = size(all_mat);

N = sum(sum(~isnan(all_mat))); % data
M = nstas+nevts;               % model parameters


si = zeros(2*N,1);
sj = zeros(2*N,1);
s = zeros(2*N,1);

d = zeros(N,1);

kk = 0;
for ie = 1:nevts
    for is = 1:nstas
        if isnan(all_mat(is,ie)), continue, end
        kk = kk+1; % n_obs
        si(2*kk+[-1 0]) = kk;
        sj(2*kk+[-1 0]) = [is, nstas+ie];
        s(2*kk+[-1 0]) = 1;
        
        d(kk) = all_mat(is,ie);
    end
end


G = sparse(si,sj,s,N,M,2*N);

L = [sparse(nevts,nstas),sparse(1:nevts,1:nevts,ones(nevts,1),nevts,nevts)];

F = [G; epsilon*L];
f = [d; zeros(nevts,1)];

m = [F'*F]\F'*f;

sta_terms = m(1:nstas);
evt_terms = m([1:nevts] + nstas);

% E = (d - G*m)'*(d - G*m)/N

% [sta_terms,nanmean(all_mat,2)]

end

