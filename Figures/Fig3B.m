%% Plot lower uptake of each metal ion
% PCA
load('sLU_res.mat');
load('CofactorYeast.mat');

fluxes = sLU_res.fluxes;
labels = sLU_res.labels;
rxns = model.rxns;
clear sLU_res;

idx_0mu = any(fluxes);
fluxes = fluxes(:,idx_0mu);
labels = labels(:,idx_0mu);
clear idx_0mu;

% remove non-metabolic fluxes
exclude = {'_translated','_cofactorbound','_enzyme_formation','_enzyme_dilution'};
metrxn = ~contains(model.rxns,exclude);
fluxes = fluxes(metrxn,:); 
rxns = rxns(metrxn,:);
clear exclude metrxn;
% remove the reactions with at least one flux > 100 or < -100
idx_over100 = any(abs(fluxes) > 100,2);
fluxes = fluxes(~idx_over100,:);
rxns = rxns(~idx_over100);
clear idx_over100;
% remove the reactions with zero flux at any simulation
idx_zero = all(fluxes == 0,2);
fluxes = fluxes(~idx_zero,:);
rxns = rxns(~idx_zero);
clear idx_zero;
% % remove the reactions with small flux (<0.0001) at any simulation
% idx_small = all(abs(fluxes) < 0.0001,2);
% fluxes = fluxes(~idx_small,:);
% rxns = rxns(~idx_small);
% clear idx_small;

% fluxes = fluxes ./ -fluxes(ismember(rxns,'r_1714'),:);


%% PCA plot
ion_id_list = {'CA';'CU';'FE';'K';'MG';'MN';'NA';'ZN'};
conserved_id_list = {'Others';'C';'N';'O';'P';'S';'Element plotted'};
conserved_color = [200,200,200% Other
                   228,26,28 % C
                   55,126,184 % N
                   77,175,74  % O
                   152,78,163  % P
                   255,127,0  % S
                   0,0,0]/255;% Key ion
alpha = [1,0.8,0.8,0.8,0.8,0.8,1];                            


labels = labels';
labels = extractBefore(labels,'_');
labels = repmat(labels,1,length(unique(ion_id_list)),1); % generate all 1 
for i = 1:length(ion_id_list)
    setother_list = setdiff(ion_id_list,ion_id_list(i));
    idx = ismember(labels(:,end),setother_list);
    labels(idx,i) = {'Others'};
end

% % tsne-2D
% Y = tsne(fluxes','Algorithm','barneshut','NumPCAComponents',50);
% 
% [~,idx2] = ismember(labels,conserved_id_list);
% idx2(idx2==0) = 7;
% colormap(conserved_color);
% 
% for i = 1:8
% subplot(2,4,i);
% gscatter(Y(:,1),Y(:,2),idx2(:,i),colormap,'.',6);legend off;
% title(ion_id_list{i},'FontSize',7,'FontName','Helvetica','Color','k');
% xlabel('tsne1','FontSize',6,'FontName','Helvetica');
% ylabel('tsne2','FontSize',6,'FontName','Helvetica');
% end

% % tsne-3D
% Y = tsne(fluxes','Algorithm','barneshut','NumPCAComponents',50,'NumDimensions',3);
% 
% [~,idx2] = ismember(labels,conserved_id_list);
% idx2(idx2==0) = 7;
% colormap(conserved_color);
% 
% for i = 1:8
% subplot(2,4,i);
% scatter3(Y(:,1),Y(:,2),Y(:,3),20,idx2(:,i),'filled');
% view(135,25);  % to make it obvious that it is a 3D plot
% title(ion_id_list{i},'FontSize',7,'FontName','Helvetica','Color','k');
% end

% %pca
figure('Name','PCA');
[~, score, ~, ~, explained, ~] = pca(fluxes','NumComponents',30);
[~,idx2] = ismember(labels,conserved_id_list);
idx2(idx2==0) = length(conserved_id_list); % key ion
for i = 1:length(ion_id_list)
    subplot(2,4,i);
    hold on;
    box on;
    for j = 1:length(conserved_id_list)
        index = find(idx2(:,i) == j);
        h = scatter(score(index,1),score(index,2),3,'o','filled','LineWidth',1,'MarkerFaceColor',conserved_color(j,:),'MarkerFaceAlpha',alpha(j));legend off;
    end
    xlim([-90 90]);
    ylim([-35 65]);
	set(gca, 'XColor','k');
    set(gca, 'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
    title(ion_id_list{i},'FontSize',7,'FontName','Helvetica','Color','k');
    xlabel(['PC1 (',num2str(round(explained(1),1)),'%)'],'FontSize',7,'FontName','Helvetica','Color','k');
    ylabel(['PC2 (',num2str(round(explained(2),1)),'%)'],'FontSize',7,'FontName','Helvetica','Color','k');
end
set(gcf,'position',[200 200 360 180]);

% colormap(conserved_color);
% for i = 1:8
% subplot(2,4,i);
% gscatter(score(:,1),score(:,2),idx2(:,i),colormap,'.',10);legend off;
% set(gca,'FontSize',6,'FontName','Helvetica');
% title(ion_id_list{i},'FontSize',7,'FontName','Helvetica','Color','k');
% xlabel(['PC1 (',num2str(explained(1)),'%)'],'FontSize',6,'FontName','Helvetica');
% ylabel(['PC2 (',num2str(explained(2)),'%)'],'FontSize',6,'FontName','Helvetica');
% end