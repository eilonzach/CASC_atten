clear all
close all
cd ~/Documents/MATLAB/CASC_atten/
addpath('matguts')

%% parms
lonlims = [-126.8 -123.3];
latlims = [ 45.5 48];
cmap = blue2red; close(gcf)
scale = 150;

keyloc = [-126.5,46];
keysiz = [0.46,0.9];

phases = {'P','S'};
components = {'Z','T'}; %'Z', 'R', or 'T'

ifsave = true;

%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% DATA DIRECTORY (top level)
datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash
% RESULTS DIRECTORY
resdir = '~/Documents/MATLAB/CASC_atten/results/'; % needs final slash


%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%
% GET EVENTS+STATIONS DATA
[ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam );
[ nstas,stas,slats,slons,selevs,ondate,offdate,staname ] = db_stadata( dbdir,dbnam );

for ip = 1:length(phases)
    phase = phases{ip};
    component = components{ip};
    tdflim = ip*1.5*[-1 1];

    % LOAD RESULTS
    load([resdir,'all_dT_',phase,'_',component,'.mat']);

    % PARSE STATIONS
    gdstas = ~isnan(nanmean(all_dT,2));
    gdstas = gdstas & slons>=lonlims(1) & slons<=lonlims(2)...
                    & slats>=latlims(1) & slats<=latlims(2);
    iobs = ~cellfun('isempty',regexp(upper(staname),'OBS '));
    itrm = ~cellfun('isempty',regexp(upper(staname),'TRAWL'));
    ilan =  cellfun('isempty',regexp(upper(staname),'OBS '));

    
    [ sta_terms,evt_terms ] = lsq_sta_evt( all_dT,0.01 );
    
    %% ------------------ MAP WITH DTSTAR FOR THIS EVENT ------------------
    figure(34), clf, hold on
    mkfig_CascMAP

    set(gcf,'position',[200 300 700 700*plot_size_ratio( lonlims,latlims )])
    axis([lonlims latlims])

    htrm = scatter(slons(gdstas&itrm),slats(gdstas&itrm),...
        scale*sqrt(sum(~isnan(all_dT(gdstas&itrm,:)),2)),sta_terms(gdstas&itrm),...
        '^','filled','MarkerEdgeColor','k');    
    hobs = scatter(slons(gdstas&iobs&~itrm),slats(gdstas&iobs&~itrm),...
        scale*sqrt(sum(~isnan(all_dT(gdstas&iobs&~itrm,:)),2)),sta_terms(gdstas&iobs&~itrm),...
        'o','filled','MarkerEdgeColor','k');
    hlan = scatter(slons(gdstas&ilan),slats(gdstas&ilan),...
        scale*sqrt(sum(~isnan(all_dT(gdstas&ilan,:)),2)),sta_terms(gdstas&ilan),...
        's','filled','MarkerEdgeColor','k');
    colormap(cmap), caxis(tdflim)

    %% colour bar
    cbar_custom(gca, 'location',[-124.5 -123.4 47.8 47.88],'tickside','bottom',...
        'lims',tdflim,'tickvals',[tdflim(1):ip*0.5:tdflim(2)],'cmap',cmap,...
        'FontSize',12,'FontWeight','bold',...
        'title',sprintf('$\\delta T_%s$ \\,(s)',phase),'interpreter','latex');

    %% key
    kle = keyloc(1) - 0.5*keysiz(1); kri = keyloc(1) + 0.5*keysiz(1);
    kbo = keyloc(2) - 0.5*keysiz(2); kto = keyloc(2) + 0.5*keysiz(2);
    kw = keysiz(1); kh = keysiz(2);

    hold on
    patch([kle kle kri kri kle],[kbo kto kto kbo kbo],...
          'w','LineWidth',2)
    
    % NEVTS
	text(keyloc(1),kto-0.03*kh,'Nevts',...
        'fontsize',15,'fontweight','bold','verticalalignment','top','horizontalalignment','center');

    scatter(kle + 0.3*kw,kbo + 0.84*kh,scale*sqrt(1),  'sk','filled','MarkerEdgeColor','k');
    scatter(kle + 0.3*kw,kbo + 0.74*kh,scale*sqrt(10), 'sk','filled','MarkerEdgeColor','k');
    scatter(kle + 0.3*kw,kbo + 0.61*kh,scale*sqrt(50),'sk','filled','MarkerEdgeColor','k');

    text(kle + 0.6*kw,kbo + 0.84*kh,'1', 'fontsize',13,'fontweight','bold','verticalalignment','middle');
    text(kle + 0.6*kw,kbo + 0.74*kh,'10','fontsize',13,'fontweight','bold','verticalalignment','middle');
    text(kle + 0.6*kw,kbo + 0.61*kh,'50','fontsize',13,'fontweight','bold','verticalalignment','middle');

    % STATYPE
	text(keyloc(1),kbo+0.52*kh,'Statype',...
        'fontsize',15,'fontweight','bold','verticalalignment','top','horizontalalignment','center');

    scatter(kle + 0.25*kw,kbo + 0.34*kh,scale*sqrt(10),  '^k','filled','MarkerEdgeColor','k');
    scatter(kle + 0.25*kw,kbo + 0.21*kh,scale*sqrt(10), 'ok','filled','MarkerEdgeColor','k');
    scatter(kle + 0.25*kw,kbo + 0.08*kh,scale*sqrt(10), 'sk','filled','MarkerEdgeColor','k');

    text(kle + 0.45*kw,kbo + 0.34*kh,'+TRM', 'fontsize',13,'fontweight','bold','verticalalignment','middle');
    text(kle + 0.45*kw,kbo + 0.21*kh,'OBS','fontsize',13,'fontweight','bold','verticalalignment','middle');
    text(kle + 0.45*kw,kbo + 0.08*kh,'Land','fontsize',13,'fontweight','bold','verticalalignment','middle');


    %% title
    title(sprintf('$\\delta T$ for $%s$-waves (%s component)',phase,component),...
          'FontSize',18,'FontWeight','bold','Interpreter','latex')
    set(gca,'FontSize',14,'LineWidth',2)
    box on
    
    % save
    if ifsave
    save2pdf(34,sprintf('FOCUS_dT_%s_%s_sta_average',phase,component),'figs');
    end
    
end % loop on phases