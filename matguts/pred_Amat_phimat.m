function [ Amat,phimat ] = pred_Amat_phimat( dtstars,dTs,A0s,fmids,alp )
% [ Amat,phimat ] = pred_Amat_Phimat( dtstars,dTs,A0s,fmids,[alpha=0] )
%   function to compute pairwise amplitude and phase spectra, given known
%   values of dtstar, dT, A0 for each station. 
if nargin < 5 || isempty(alp)
    alp = 0;
end

Nstas = length(dtstars);
Npair = handshake(Nstas);
Nf = length(fmids);

Amat = zeros(Npair,Nf);
phimat = zeros(Npair,Nf);

ws = 2*pi*fmids;

count = 0;
for is1 = 1:Nstas
for is2 = is1+1:Nstas
    count = count+1;
    
    RA0 = A0s(is2)/A0s(is1);
    dtst = dtstars(is2)-dtstars(is1);
    dT = dTs(is2)-dTs(is1);
    
    if alp == 0;
        lnApred   = log(RA0) - pi*fmids*dtst;
        phipred = -(1/pi)*log(fmids)*dtst + dtst + dT;
    else
        lnApred = log(RA0) - 0.5*ws.^(1-alp)*dtst;
        phipred = 0.5*ws.^(-alp).*(2*pi).^alp * cot(alp*pi/2)*dtst + dT;
    end
    
    Amat(count,:) = exp(lnApred);
    phimat(count,:) = phipred;    

end
end



end

