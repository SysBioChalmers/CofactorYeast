%% Plot iron

load('sI_res.mat');
load('CofactorYeast.mat');
load('enzymedata.mat');
load('CofactorDataset.mat');

selected_data = 1:2:19;
fluxes = [sI_res.fluxes(:,selected_data) sI_res.flux_ref];
label_tmp = num2str(sI_res.k_cf(1,selected_data));
label_tmp = strsplit(label_tmp);
label = [label_tmp 'Ref'];

%% Exchange fluxes
glc = -1*fluxes(strcmp(model.rxns,'r_1714'),:);
etoh = fluxes(strcmp(model.rxns,'r_1761'),:);
mu = fluxes(strcmp(model.rxns,'r_2111'),:);
fe = -1*fluxes(strcmp(model.rxns,'r_1861'),:);
glyc = fluxes(strcmp(model.rxns,'r_1808'),:);

figure('Name','1');
subplot(4,1,1);
b1 = bar(mu,0.7,'FaceColor','flat','LineWidth',0.1);
b1.CData(1:length(selected_data),:) = repmat([242,94,13]/255,length(selected_data),1);
b1.CData(end,:) = [0 0 0];
b1.EdgeColor = 'w';
xticks(1:1:length(label)+1);
xticklabels(label);
xtickangle(90);
ylim([0.28 0.42]);
% yticks(0.3:0.02:0.48);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica','Color','k');
subplot(4,1,2);
b2 = bar(glc,0.7,'FaceColor','flat','LineWidth',0.1);
b2.CData(1:length(selected_data),:) = repmat([242,94,13]/255,length(selected_data),1);
b2.CData(end,:) = [0 0 0];
b2.EdgeColor = 'w';
xticks(1:1:length(label)+1);
xticklabels(label);
xtickangle(90);
ylim([6 24]);
set(gca,'FontSize',6,'FontName','Helvetica');
title('Glucose uptake','FontSize',7,'FontName','Helvetica','Color','k');
subplot(4,1,3);
b3 = bar(etoh,0.7,'FaceColor','flat','LineWidth',0.1);
b3.CData(1:length(selected_data),:) = repmat([242,94,13]/255,length(selected_data),1);
b3.CData(end,:) = [0 0 0];
b3.EdgeColor = 'w';
xticks(1:1:length(label)+1);
xticklabels(label);
xtickangle(90);
ylim([16 34]);
set(gca,'FontSize',6,'FontName','Helvetica');
title('Ethanol production','FontSize',7,'FontName','Helvetica','Color','k');
ylabel('Flux (mmol/gCDW/h)','FontSize',7,'FontName','Helvetica','Color','k');
subplot(4,1,4);
b4 = bar(glyc,0.7,'FaceColor','flat','LineWidth',0.1);
b4.CData(1:length(selected_data),:) = repmat([242,94,13]/255,length(selected_data),1);
b4.CData(end,:) = [0 0 0];
b4.EdgeColor = 'w';
xticks(1:1:length(label)+1);
xticklabels(label);
xtickangle(90);
ylim([0 5]);
set(gca,'FontSize',6,'FontName','Helvetica');
title('Glycerol production','FontSize',7,'FontName','Helvetica','Color','k');
xlabel('theta','FontSize',7,'FontName','Helvetica','Color','k');

set(gcf,'position',[300 300 120 300]);

%% Protein levels
protein_conc_ref = calculateProteinConc(model,model.genes,fluxes(:,end));
protein_conc_low = calculateProteinConc(model,model.genes,fluxes(:,1:end-1);
protein_conc = [protein_conc_low protein_conc_ref];

% remove low absolute protein level in reference
cutoff_low_abs = 0.05;
low_abs_value = quantile(protein_conc_ref(protein_conc_ref>0),cutoff_low_abs);
protein_list = model.genes(protein_conc_ref>low_abs_value);
data_abs = protein_conc(protein_conc_ref>low_abs_value,:);

data_rel = data_abs(:,1:end-1)./data_abs(:,end);

[~,txt,~] = xlsread('Yeast8_Modification.xlsx','SGDgeneNames');
gname_1 = txt(2:end,1);
gname_2 = txt(2:end,2);
clear txt;
gname_2(ismember(gname_2,'')) = gname_1(ismember(gname_2,''));

[~,b] = ismember(protein_list,gname_1);
protein_id_list = gname_2(b);









