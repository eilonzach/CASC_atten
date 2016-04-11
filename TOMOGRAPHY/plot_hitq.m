function plot_hitq( par )
% plot_hitq( par )
%   plot the hit quality

shape = [par.ny,par.nx,par.nz];

lats= reshape(par.mlt,shape);
lons= reshape(par.mln,shape);
hitq = reshape(par.hitq,shape);
hitc = reshape(sum(par.hitc,2),shape);

figure(66), clf
for iz = 2:par.nz
    subplot(5,3,iz-1)
    contourf(lons(:,:,iz),lats(:,:,iz),100*hitq(:,:,iz),[0:10:100])
%     contourf(lons(:,:,iz),lats(:,:,iz),100*hitc(:,:,iz),[0:10:100])
%     xlabel('Longitude'); 
%     ylabel('Latitude');
    title(sprintf('Depth %.0f km',par.zz(iz)),'FontSize',14)
    xlim([-132 -120]); ylim([39 51])
    shading flat

    set(gcf,'Position',[0,0,800,800])
end

end

