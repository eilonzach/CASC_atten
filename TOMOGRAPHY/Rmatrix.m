function [ r_k ] = Rmatrix( F, par, pt_k )
%[ r_k ] = Rmatrix( F, par, pt_k )
%   compute a row of the resolution matrix corresponding to model parameter
%   element k - chosen as the model parameter nearest to pt_k in 3D space
%   i.e. pt_k = [lt ln z];

k = mindex(abs(par.mlt-pt_k(1)) + abs(par.mln-pt_k(2)) + abs(par.mz-pt_k(3)));

% R = F'*F;
m_k = zeros(size(F,2),1);
m_k(k) = 1;

f = F*m_k;

[r_k,~,~,~,~] = lsqr( F, f, 1e-4, 400 );
% r_k = full(diag(F'*F));
% r_k = full(diag(G'*G));
% [r_k,~,~,~,~] = lsqr( G, G*m_k, 1e-4, 400 );
% [r_k,~,~,~,~] = lsqr( H_damp, H_damp*m_k, 1e-4, 400 );
% [r_k,~,~,~,~] = lsqr( H_smth, H_smth*m_k(1:par.nmodel), 1e-4, 400 );
figure(8);
plot(r_k)

model.mval = r_k(1:par.nmodel)*1e2;
[plot_rk] = conv2plotable(model,par);
plot_basic(plot_rk,par,1,0);


end

