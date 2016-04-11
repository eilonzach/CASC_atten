function [ synth_model ] = make_synth_model( par )

if strcmp(par.sym.opt,'custom')

    mvv = par.mvav;
    mva = zeros(size(mvv));

    na =  length(par.sym.aa);
    for ia = 1:na
        ind1 = par.mx <= par.sym.acx(ia) + 0.5*par.sym.awx(ia);
        ind2 = par.mx >= par.sym.acx(ia) - 0.5*par.sym.awx(ia);
        ind3 = par.my <= par.sym.acy(ia) + 0.5*par.sym.awy(ia);
        ind4 = par.my >= par.sym.acy(ia) - 0.5*par.sym.awy(ia);
        ind5 = par.mz <= par.sym.acz(ia) + 0.5*par.sym.awz(ia);
        ind6 = par.mz >= par.sym.acz(ia) - 0.5*par.sym.awz(ia);
        ind = ind1 & ind2 & ind3 & ind4 & ind5 & ind6;

        mvv(ind) = mvv(ind)*(1 + 0.01*par.sym.adv(ia));
        mva(ind) = mva(ind) + par.sym.aa(ia);
    end

    [mvx,mvy] = va2xy(mvv,mva,'forward');

    %add a little random noise
    mvx = mvx + random('norm',0,0.00001,size(mvx));
    mvy = mvy + random('norm',0,0.00001,size(mvy));

    synth_model.mvx = mvx;
    synth_model.mvy = mvy;
    
elseif strcmp(par.sym.opt,'checker')

    synth_model = make_checker_model(par,par.sym.dq);

end



end

