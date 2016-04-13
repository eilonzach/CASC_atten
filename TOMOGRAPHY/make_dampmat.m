function [ H_damp,h_damp ] = make_dampmat( data,par )
% makes damping matrix 
% push layer average to zero
% push static terms to zero
% 
% model parameter vector is [nmodel,nevts,nstas]

fprintf('creating damping matrix... \n');

damp_lav0 = 10; % damping layer average to zero

nstas = data.stn.nstas;   
nevts = data.evt.nevts;
nmod = par.nmodel;                          
nmat = nmod+nevts+nstas; % horizontal dimension of smoothing matrix

zz = par.zz;

%% OVERALL damping (i.e. x,y=0)
% make weights for diagonal - damping all to zeros
mwt = ones(nmod,1);

% ewt = par.damp_evt*ones(2*nevts,1);
% swt = par.damp_stn*ones(2*nstas,1);
ewt = ones(nevts,1);
swt = ones(nstas,1);

wt = [mwt;  
      ewt; 
      swt];    % repeat mwt for x,y at each node

% turn into sparse diagonal matrix
si = 1:nmat;
sj = 1:nmat;
s  = 0;
fprintf('NB - NO DIAGONAL DAMPING\n')

H_damp = sparse(si,sj,s,nmat,nmat);
h_damp = zeros(nmat,1);

%% add event damping 
mwt = ones(nevts,1);
eadi = [1:nevts]';
eadj = [1:nevts]'+nmod;
ead  = [mwt];
ED = sparse(eadi,eadj,ead,nevts,nmat);

H_damp = [H_damp;ED];
h_damp = [h_damp;par.start_model.estatic];


%% add station damping 
mwt = ones(nstas,1);
sadi = [(1:nstas)]';
sadj = [(1:nstas)]'+nmod+nevts;
sad  = [mwt];
SD = sparse(sadi,sadj,sad,nstas,nmat);

H_damp = [H_damp;SD];
h_damp = [h_damp;par.start_model.sstatic];

%% add constraint that average in each layer is 0
bi = zeros(nmod,1);
bj = zeros(nmod,1);
b  = zeros(nmod,1);
nlay = par.nx*par.ny;
n = 1;
for iz = 1:par.nz
    bi(n:n+nlay-1) = iz*ones(nlay,1);
    bj(n:n+nlay-1) = find(par.mz==zz(iz));
    b(n:n+nlay-1)  = ones(nlay,1);
    n = n+nlay;
end
B = sparse(bi,bj,b,par.nz,nmat);

H_damp = [H_damp;damp_lav0*B];
h_damp = [h_damp;zeros(par.nz,1)];
    
%% Smooth by ocean floor age
% smooth in age bins of 1 Ma. Only do for ages of < 12 Ma
if par.age_smooth
disp('Smoothing along isochrons ');
    
ages = 0:1:12;
Nag = length(ages)-1;
ci = nan(2*nmod,1); % just to be safe, max # of elements is 2*nmod
cj = nan(2*nmod,1);
c  = nan(2*nmod,1);

n = 1; % number of elements in sparse index vectors
iazi = 0; % number of age-depth perms
for iz = 1:par.nz
for ia = 1:Nag
    % find points in this depth slice within the age bounds
    inds = find(par.mage <= ages(ia+1) & par.mage >= ages(ia) & par.mz == zz(iz)); 
    ninds = length(inds);
    inds = [inds(:);inds(1)]; % just tack on first element at the end to complete loop

    %     Now string these all together by looping through them minimising
    %     first difference between each and next
    for ii = 1:ninds
        iazi = iazi+1;
        
        ci(n+[0:1]) = iazi*[1 1];
        cj(n+[0:1]) = [inds(ii) inds(ii+1)];
        c(n+[0:1])  = [-1 1] ;
        
        n = n+2;
    end
end
end
ci(isnan(ci))=[];
cj(isnan(cj))=[];
c(isnan(c))=[];

C = sparse(ci,cj,c,iazi,nmat);

H_damp = [H_damp;10*C];
h_damp = [h_damp;zeros(iazi,1)];

end

% %% damp out any anisotropy in the south of the array 
% % since splitting gets complex here
% if par.no_S_anis
%     
% indx = find(par.mlt < -10);
% nindx = length(indx);
% di = [(1:nindx),(1:nindx)]';
% dj = [indx,nmod + indx];
% d  = [ones(nindx,1);-ones(nindx,1)];
% 
% D = sparse(di,dj,d,nindx,nmat);
% 
% H_damp = [H_damp;1e4*D];
% 
% end

end

    

