%% Plot Biolog
load('sB_res.mat');

load('CofactorYeast.mat');
load('enzymedata.mat');

ion_id_list = {'CA';'CU';'FE';'K';'MG';'MN';'NA';'ZN'};
ion_ex_list = {'r_4600'; ... % Ca(2+) exchange
               'r_4594'; ... % Cu(2+) exchange
               'r_1861'; ... % iron(2+) exchange
               'r_2020'; ... % potassium exchange
               'r_4597'; ... % Mg(2+) exchange
               'r_4595'; ... % Mn(2+) exchange
               'r_2049'; ... % sodium exchange
               'r_4596'};... % Zn(2+) exchange

%% Predicted ion usage
[~,txt,~] = xlsread('Yeast8_Modification.xlsx','SGDgeneNames');
gname_1 = txt(2:end,1);
gname_2 = txt(2:end,2);
clear txt;
gname_2(ismember(gname_2,'')) = gname_1(ismember(gname_2,''));

[~,b] = ismember(model.genes,gname_1);
protein_id_list = gname_2(b);
clear b;

gname_2 = lower(gname_2);
gname_2 = cellfun(@(x) strcat(upper(x(1)),x(2:end)),gname_2,'UniformOutput',false);

top = 5;

figure('Name','1');
for i = 1:length(ion_id_list)
    ionusage_tmp = calculateCofactorUsage4protein(model,ion_id_list{i},model.genes,sB_res.fluxes);
    idx_ion_ex_tmp = ismember(model.rxns,ion_ex_list(i));
    ion_uptake_tmp = -1*sB_res.fluxes(idx_ion_ex_tmp,:);
    tot_tmp = ion_uptake_tmp ./ sB_res.fluxes(strcmp(model.rxns,'r_2111'),:);
    perc_list_tmp = ionusage_tmp./tot_tmp*100;
    median_tmp = median(perc_list_tmp,2);
    
    
    [~,I] = sort(median_tmp,'descend');
    proteinid = model.genes(I(1:top));
    value = perc_list_tmp(I(1:top),:);
    [~,b] = ismember(proteinid,gname_1);
    proteinname = gname_2(b);
    proteinname(median(value,2) == 0) = {'n/a'};
    
    subplot(2,4,i);
    hold on;
    h = boxplot(value','Symbol','.','OutlierSize',8,'Widths',0.5,'Colors',[242,94,13]/255);
    set(h,{'linew'},{0.5});
    xlim([0.1 length(proteinname)+0.9]);
    xticks(1:1:length(proteinname));
    set(gca, 'XColor','k');
    set(gca, 'YColor','k');
    set(gca,'XTickLabel',proteinname);
    set(gca,'FontSize',6,'FontName','Helvetica');
	set(gca,'XColor','k');
    set(gca,'YColor','k');
    ylabel('Usage fraction (%)','FontSize',7,'FontName','Helvetica','Color','k');
    
    xtickangle(90);
    title(ion_id_list{i},'FontSize',7,'FontName','Helvetica','Color','k');
    box off;
end
set(gcf,'position',[200 500 350 180]);




