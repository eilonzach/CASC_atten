function [ checker_model ] = make_checker_model(par,dval)

mf = zeros(par.ny,par.nx,par.nz);


%% horiz spacing for upper layers
mfx = [ 0 0 1 1 1 1 0 0 -1 -1 -1 -1];
mfx = repmat(mfx,1,30);
mfx1 = mfx(1:par.nx);

mfy = [ 0 0 1 1 1 1 0 0 -1 -1 -1 -1]';
mfy = repmat(mfy,30,1);
mfy1 = mfy(1:par.ny);

mfxy1 = mfy1*mfx1; % 2D grid

%% horiz spacing for lower layers
mfx = [ 0 0 1 1 1 1 0 0 -1 -1 -1 -1];
mfx = repmat(mfx,1,30);
mfx2 = mfx(1:par.nx);

mfy = [ 0 0 1 1 1 1 0 0 -1 -1 -1 -1]';
mfy = repmat(mfy,30,1);
mfy2 = mfy(1:par.ny);

mfxy2 = mfy2*mfx2; % 2D grid

%% vertical spacing
mfz = [0 0  1  1  0 -1 -1  0  0 ];
mfz = repmat(mfz,1,2);
mfz = mfz(1:par.nz);

%% now rep over layers
for iz = 1:par.nz
    if iz < 5
        mfxy = mfxy1; 
    else
        mfxy = mfxy2;
    end
    mf(:,:,iz) = mfxy*mfz(iz);
end

mval = (0.01*dval)*mf(:); % 0.01 to convert from percent


%add a little random noise
mval = mval + random('norm',0,0.00001,size(mval));

checker_model.mval = mval;
   
