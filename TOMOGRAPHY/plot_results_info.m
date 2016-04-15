function plot_results_info(par,model,data,d_use,res,start_model)
% plot some of the latest model information

M = par.nmodel;
nevts = data.evt.nevts;
nstas = data.stn.nstas;


figure(57), clf
set(gcf,'Position',[10  10   930   780]);

% plot event terms
subplot(14,5,[1:3,6:8,11:13,16:18])
te = model.estatic; 
[Ne,xe]=hist(-te); % NEGATIVE, so that positive value on plot means more time beneath model
bar(xe,Ne,'stack');
xlim([-3 3])
xlabel('Event static times','FontSize',14)
ylabel('Number','FontSize',14)
% title('Evt static','FontSize',16)

% plot station terms
subplot(14,5,[26:28,31:33,36:38,41:43])
tc = model.sstatic;
[Nc,xc]=hist(-tc); % NEGATIVE, so that positive value on plot means more time beneath model
bar(xc,Nc,'stack');
xlim([-3 3])
xlabel('Station static times','FontSize',14)
ylabel('Number','FontSize',14)
% title('Stn static','FontSize',16)

% plot resid vs. dT for differential travel times
subplot(14,5,[51:53,56:58,61:63,66:68])
hold on
plot([-3 3],0*[1 1],'--r','LineWidth',2)
scatter(d_use,res,10*data.ray.wt)
axis([-2.5 2.5 -2 2])
xlabel('Data (s)','FontSize',14)
ylabel('Residual (s)','FontSize',14)
% title('Resid vs. data ','FontSize',16)

% locate event terms
subplot(14,5,[4,5,9,10,14,15,19,20,24,25,29,30])
for ie = 1:data.evt.nevts
    rp(ie) = mean(data.ray.pd(data.ray.orid==data.evt.orid(ie)));
    bz(ie) = mean(data.ray.baz(data.ray.orid==data.evt.orid(ie)));
end
scatter(bz,rp,70,model.estatic,'filled')
xlim([0 360])
caxis([-3 3]);
xlabel('Event static times','FontSize',14)
% ylabel('Number','FontSize',14)
% title('Stn static','FontSize',16)

% map station terms
subplot(14,5,[39,40,44,45,49,50,54,55,59,60,64,65,69,70])
scatter(data.stn.lon,data.stn.lat,70,model.sstatic,'filled')
axis([-131 -120 39 51])
xlabel('Station static times','FontSize',14)
caxis([-2 2])
% ylabel('Number','FontSize',14)
% title('Stn static','FontSize',16)



if nargin == 6 && ~isempty(start_model)
    figure(58), clf
    set(gcf,'Position',[640  10   630   830]);

    % plot event terms
    subplot(3,1,1)
    dte = model.estatic-start_model.estatic;
    [Ne,xe]=hist(dte); % NEGATIVE, so that positive value on plot means more time beneath model
    bar(xe,Ne,'stack');
    xlim([-3 3])
    xlabel('Error (s)','FontSize',14)
    ylabel('Number','FontSize',14)
    title('Evt static errors','FontSize',16)

    % plot station terms
    subplot(3,1,2)
    dtc = model.sstatic-start_model.sstatic;
    [Nc,xc]=hist(dtc); % NEGATIVE, so that positive value on plot means more time beneath model
    bar(xc,Nc,'stack');
    xlim([-3 3])
    xlabel('Error (s)','FontSize',14)
    ylabel('Number','FontSize',14)
    title('Station static errors','FontSize',16)

    % plot resid vs. dT for differential travel times
    subplot(3,1,3)
    hold on
    plot(0.05*[1 1],0.05*[1 1],'--r','LineWidth',2)
    plot(start_model.mval,model.mval,'o')
    axis([-0.06 0.06 -0.06 0.06])
    xlabel('Start model','FontSize',14)
    ylabel('Out model','FontSize',14)
    title('In vs. out model ','FontSize',16)
    
end
return


%% ==============    more detail on event and station terms   ==============    
figure(58)
set(gcf,'Position',[100  10   1200   500]);
clf

% Plot station terms at their lat, lon
subplot(1,2,1), hold on
% NEGATIVE SIGNS BELOW, so that positive value on plot means more time in crust
scatter(data.stn.lon,data.stn.lat,65*ones(nstas,1),-tc,'fill')
scatter(data.stn.lon,data.stn.lat,150*ones(nstas,1),-dtc,'Linewidth',2)% this value is positive if tEcrust > tNcrust
% plot(coast.ncst,'LineWidth',2)
xlabel('Lon','FontSize',14)
ylabel('Lat','FontSize',14)
title('Stn static','FontSize',16)

hc2 = colorbar; caxis([-1.5 1.5])
set(get(hc2,'Ylabel'),'string','Centre=Tav, Edge=dt','Fontsize',16,'FontWeight','bold');


% Plot event terms by inc and seaz
subplot(1,2,2)
for ie = 1:data.evt.nevts
    ind = data.ray.orid == data.evt.orid(ie);
    rayp(ie) = mean(data.ray.pd(ind));
    gcarc(ie) = mean(data.ray.gcarc(ind));
    seaz(ie) = mean(data.ray.baz(ind));
    nars(ie) = double(sum(ind(ind)./data.ray.sd(ind).^2));
end
h = scatter(seaz,rayp,0.1*nars,-dte,'LineWidth',2);

xlim([00 360])
ylim([10 16])
xlabel('seaz','FontSize',14), set(gca,'XTick',[0:45:360]);
ylabel('rayp (s/d)','FontSize',14)
title('Evt static','FontSize',16)

hc1 = colorbar; caxis([-1.5 1.5])
set(get(hc1,'Ylabel'),'string','<<< N-fast     Tn-Te     E-fast >>>','Fontsize',16,'FontWeight','bold');



end