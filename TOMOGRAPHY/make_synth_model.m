function [ synth_model ] = make_synth_model( par )

if strcmp(par.sym.opt,'custom')

    mval = zeros(size(par.nmodel));

    na =  length(par.sym.aa);
    for ia = 1:na
        ind1 = par.mx <= par.sym.acx(ia) + 0.5*par.sym.awx(ia);
        ind2 = par.mx >= par.sym.acx(ia) - 0.5*par.sym.awx(ia);
        ind3 = par.my <= par.sym.acy(ia) + 0.5*par.sym.awy(ia);
        ind4 = par.my >= par.sym.acy(ia) - 0.5*par.sym.awy(ia);
        ind5 = par.mz <= par.sym.acz(ia) + 0.5*par.sym.awz(ia);
        ind6 = par.mz >= par.sym.acz(ia) - 0.5*par.sym.awz(ia);
        ind = ind1 & ind2 & ind3 & ind4 & ind5 & ind6;

        mval(ind) = mval(ind)*(1 + 0.01*par.sym.adval(ia));
    end

    %add a little random noise
    mval = mval + random('norm',0,0.00001,size(mval));

    synth_model.mval = mval;
    
elseif strcmp(par.sym.opt,'checker')

    synth_model = make_checker_model(par,par.sym.dval);

end



end

