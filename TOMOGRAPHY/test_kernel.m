% pause
kval = zeros(size(K.n_indx));
klt = zeros(size(K.n_indx));
kln = zeros(size(K.n_indx));
for ik = 1:length(K.n_vals)
   kval(ik) = sum(K.n_vals{ik});
   klt(ik) = mean(par.mlt(K.n_indx{ik}));
   kln(ik) = mean(par.mln(K.n_indx{ik}));
   if kval(ik)<20, return, end
end

figure(1), clf,hold on
scatter(kln,klt,40,kval,'filled')
colorbar
plot(data.stn.lon,data.stn.lat,'v')


    
return


kval = zeros(size(par.mz));
figure(1), clf, hold on
gd = par.mln<-120 & par.mln>-132 & par.mlt>38 & par.mlt<51;
for ik = 1:length(K.n_vals)
    ind = false(par.nmodel,1);
    ind(K.n_indx{ik}) = true;
    kval(ind) = kval(ind) + K.n_vals{ik};
    set(gca,'zdir','reverse','xlim',[-132 -120],'ylim',[39 50])
    clf, hold on
    scatter3(par.mln(gd),par.mlt(gd),par.mz(gd),20,kval(gd),'filled')
    h = scatter3(par.mln(gd&ind),par.mlt(gd&ind),par.mz(gd&ind),120,kval(gd&ind));
    pause
end
    