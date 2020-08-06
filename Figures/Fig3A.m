%% Plot lower uptake of each metal ion
% plots
load('sLU_res.mat');
load('CofactorYeast.mat');

%%
ion_id_list = {'CA';'CU';'FE';'K';'MG';'MN';'NA';'ZN'};

figure('Name','plot');
for i = 1:length(ion_id_list)
    ion = ion_id_list{i};
    idx = contains(sLU_res.labels,ion);
    fluxes_tmp = sLU_res.fluxes(:,idx);
    labels_tmp = sLU_res.labels(1,idx);
    labels_tmp = cellfun(@(x) x(strfind(x,'_')+1:end),labels_tmp,'UniformOutput',false);
    labels_tmp = strrep(labels_tmp,'_','.');
    labels_tmp = cellfun(@(x) str2double(x),labels_tmp,'UniformOutput',false);
    lower_values = cell2mat(labels_tmp);
    
    display([num2str(i),'/',num2str(length(ion_id_list))]);
    
    mu = fluxes_tmp(ismember(model.rxns,'r_2111'),:);
    glc = -1*fluxes_tmp(ismember(model.rxns,'r_1714'),:);
    
    subplot(2,4,i);
    
    hold on;
    plot(lower_values*100,mu,'-o','MarkerSize',0.5,'LineWidth',0.75,'Color',[242,94,13]/255);
    set(gca,'FontSize',6,'FontName','Helvetica');
    set(gca,'ycolor','k');
    ylim([0 0.45]);
    
    if i == 1 || i == 5
        ylabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica','Color','k');
    end
    
    title(ion_id_list{i},'FontSize',7,'FontName','Helvetica','Color','k');
    xlabel('Relative uptake (%)','FontSize',7,'FontName','Helvetica','Color','k');
end
set(gcf,'position',[200 500 280 100]);



% for i = 1:length(ion_id_list)
%     ion = ion_id_list{i};
%     idx = contains(sLU_res.labels,ion);
%     fluxes_tmp = sLU_res.fluxes(:,idx);
%     labels_tmp = sLU_res.labels(1,idx);
%     labels_tmp = cellfun(@(x) x(strfind(x,'_')+1:end),labels_tmp,'UniformOutput',false);
%     labels_tmp = strrep(labels_tmp,'_','.');
%     labels_tmp = cellfun(@(x) str2double(x),labels_tmp,'UniformOutput',false);
%     lower_values = cell2mat(labels_tmp);
%     
%     display([num2str(i),'/',num2str(length(ion_id_list))]);
%     
%     mu = fluxes_tmp(ismember(model.rxns,'r_2111'),:);
%     glc = -1*fluxes_tmp(ismember(model.rxns,'r_1714'),:);
%     etoh = fluxes_tmp(ismember(model.rxns,'r_1761'),:);
%     yield = mu./(glc*180/1000);
%     
%     subplot(2,4,i);
%     
%     hold on;
%     yyaxis left;
%     plot(lower_values,mu,'-o','MarkerSize',1,'LineWidth',0.75,'Color','k');
%     set(gca,'FontSize',6,'FontName','Helvetica');
%     set(gca,'ycolor','k');
%     ylim([0 0.45]);
%     
%     if i == 1 || i == 5
%         ylabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica','Color','k');
%     end
%     
%     yyaxis right;
%     plot(lower_values,yield,'-o','MarkerSize',1,'LineWidth',0.75,'Color',[242,94,13]/255);
%     set(gca,'ycolor','r');
%     ylim([0 0.45]);
%     
%     if i == 4 || i == 8
%         ylabel('Biomass yield','FontSize',7,'FontName','Helvetica','Color',[242,94,13]/255);
%     end
%     
%     title(ion_id_list{i},'FontSize',7,'FontName','Helvetica','Color','k');
%     xlabel('Relative uptake','FontSize',7,'FontName','Helvetica','Color','k');
% end
% set(gcf,'position',[200 500 300 130]);


