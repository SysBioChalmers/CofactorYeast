%% Plot lower uptake
load('sLU_res.mat');

load('CofactorYeast.mat');
load('enzymedata.mat');

%% Data
lbl = sLU_res.labels;
lbl = cellfun(@(x) x(1:strfind(x,'_')-1),lbl,'UniformOutput',false);
lbl = unique(lbl);

figure('Name','1');

for i = 1:length(lbl)
    ion = lbl{i};
    idx = contains(sLU_res.labels,ion);
    fluxes_tmp = sLU_res.fluxes(:,idx);
    labels_tmp = sLU_res.labels(1,idx);
    labels_tmp = cellfun(@(x) x(strfind(x,'_')+1:end),labels_tmp,'UniformOutput',false);
    labels_tmp = strrep(labels_tmp,'_','.');
    labels_tmp = cellfun(@(x) str2double(x),labels_tmp,'UniformOutput',false);
    lower_values = cell2mat(labels_tmp);
    mu = fluxes_tmp(ismember(model.rxns,'r_2111'),:);
    glc = -1*fluxes_tmp(ismember(model.rxns,'r_1714'),:);
    etoh = fluxes_tmp(ismember(model.rxns,'r_1761'),:);
    yield = mu./(glc*180/1000);
    
    plotidx = 1:1:length(lower_values);
    
    subplot(3,length(lbl),i);
    hold on;
    box on;
    plot(lower_values(plotidx),mu(plotidx),'-o','MarkerSize',1,'LineWidth',0.75,'Color',[140,81,10]/255);
    xlim([0 1]);
    ylim([0 0.45]);
    set(gca, 'XColor','k');
    set(gca, 'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
    if i == 1
        ylabel('Growth (/h)','FontSize',6,'FontName','Helvetica');
    end
    title(ion,'FontSize',7,'FontName','Helvetica','Color','k');
    
    subplot(3,length(lbl),i+length(lbl));
    hold on;
    box on;
    plot(lower_values(plotidx),glc(plotidx),'-o','MarkerSize',1,'LineWidth',0.75,'Color',[55,126,184]/255);
    plot(lower_values(plotidx),etoh(plotidx),'-o','MarkerSize',1,'LineWidth',0.75,'Color',[255,127,0]/255);
    xlim([0 1]);
    ylim([0 40]);
    set(gca, 'XColor','k');
    set(gca, 'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
    if i == 1
        ylabel('Exchange rate','FontSize',6,'FontName','Helvetica');
    end
    
    subplot(3,length(lbl),i+2*length(lbl));
    hold on;
    box on;
    plot(lower_values(plotidx),yield(plotidx),'-o','MarkerSize',1,'LineWidth',0.75,'Color',[1,102,94]/255);
    xlim([0 1]);
    ylim([0 0.55]);
    set(gca, 'XColor','k');
    set(gca, 'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
    if i == 1
        ylabel('Biomass yield','FontSize',6,'FontName','Helvetica');
    end
    xlabel('Relative uptake','FontSize',6,'FontName','Helvetica');
    
end

set(gcf,'position',[300 500 550 180]);







