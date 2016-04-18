function plot_hitq( par,opt )
% plot_hitq( par, [opt=1])
%   plot the hit quality (opt=1) or count (opt=2)


if nargin < 2 || isempty(opt)
    opt=1;
end

shape = [par.ny,par.nx,par.nz];

lats= reshape(par.mlt,shape);
lons= reshape(par.mln,shape);
hitq = reshape(par.hitq,shape);
hitc = reshape(sum(par.hitc,2),shape);

if opt == 1
    hit = 100*hitq;
    str = 'Hit Quality ';
elseif opt == 2
    hit = hitc;
    str = 'Hit Count ';
end
datstr= {'$\delta$ T ','$\Delta \, t^*$ '};    
compstr= {'P','S'};    

figure(66), clf
    set(gcf,'Position',[0,0,800,800])
for iz = 2:par.nz
    subplot(4,3,iz-1)
    contourf(lons(:,:,iz),lats(:,:,iz),hit(:,:,iz),[0:10:100])
%     xlabel('Longitude'); 
%     ylabel('Latitude');
    title(sprintf('Depth %.0f km',par.zz(iz)),'FontSize',14)
    xlim([-132 -120]); ylim([39 51])
    shading flat 
    caxis([0 100])
    
    
    % label the type of data
    if iz==par.nz-1
        xlabel([str,datstr{par.t_ts},compstr{par.PS}],...
            'interpreter','latex','FontSize',20)
    end
    

end

end

