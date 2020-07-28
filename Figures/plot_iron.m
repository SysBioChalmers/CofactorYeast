load('sI_res.mat');
load('CofactorYeast.mat');
load('enzymedata.mat');
load('CofactorDataset.mat');

%% plot fluxes
fluxes = [sI_res.fluxes sI_res.flux_ref];
glc = -1*fluxes(strcmp(model.rxns,'r_1714'),:);
etoh = fluxes(strcmp(model.rxns,'r_1761'),:);
mu = fluxes(strcmp(model.rxns,'r_2111'),:);
fe = -1*fluxes(strcmp(model.rxns,'r_1861'),:);


figure('Name','1');
label_tmp = num2str(sI_res.k_cf);
label_tmp = strsplit(label_tmp);
label = [label_tmp 'Ref'];
subplot(3,1,1);
b1 = bar(glc,'FaceColor','flat','LineWidth',1);
b1.CData(end,:) = [1 1 1];
xticks(1:1:length(label)+1);
xticklabels(label);
ylim([12 21]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Glucose (mmol/gCDW/h)','FontSize',7,'FontName','Helvetica');
subplot(3,1,2);
b2 = bar(etoh,'FaceColor','flat','LineWidth',1);
b2.CData(end,:) = [1 1 1];
xticks(1:1:length(label)+1);
xticklabels(label);
ylim([20 36]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Ethanol (mmol/gCDW/h)','FontSize',7,'FontName','Helvetica');
subplot(3,1,3);
b3 = bar(mu,'FaceColor','flat','LineWidth',1);
b3.CData(end,:) = [1 1 1];
xticks(1:1:length(label)+1);
xticklabels(label);
ylim([0.3 0.48]);
yticks(0.3:0.02:0.48);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[300 300 130 330]);


%% plot protein changes

protein_conc_ref = calculateProteinConc(model,model.genes,sI_res.flux_ref);
protein_conc_low = calculateProteinConc(model,model.genes,sI_res.fluxes);
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

% figure('Name','2');
% data = log2(data_rel);
% maxvalue = max(max(data));
% minvalue = min(min(data));
% data(data == -inf) = -maxvalue;
% 
% h = heatmap(label_tmp,protein_id_list,data,'Colormap',jet,'ColorMethod','count','CellLabelColor','none');
% h.Title = 'log2FC';
% set(h,'FontSize',6,'FontName','Helvetica');
% set(gcf,'position',[500 300 700 120]);
% set(gca,'position',[0.03 0.4 0.9 0.5]);

figure('Name','3');
color_set = [153,153,153
             197,86,89
             91,183,205
             203,180,123
             117,184,181
             71,120,185
             84,172,117
             220,220,220]/255;
top_proteins = 6;
plotidx = 1:1:length(sI_res.k_cf);
load('cofactor_info.mat');

cfusage = calculateCofactorUsage4protein(model,'FE',model.genes,fluxes);

tot_tmp = fe ./ mu;
perc_list_tmp = cfusage./tot_tmp*100;

idx_ion = ismember(cofactor_info.element_id,'FE');
tot_unmodeled_tmp = cofactor_info.element_abund_total(idx_ion) - cofactor_info.element_abund_modeled(idx_ion);
tot_unmodeled_tmp = tot_unmodeled_tmp * 1e3 / (6.02e23*13e-12);

max_tmp = max(perc_list_tmp(:,1),[],2);

if length(find(max_tmp)) >= top_proteins
    top = top_proteins;
else
    top = length(find(perc_list_tmp(:,1)));
end

[~,I] = sort(max_tmp,'descend');
proteinid = model.genes(I(1:top));
value = perc_list_tmp(I(1:top),:);
perc_unmodeled_fe = tot_unmodeled_tmp./tot_tmp*100;
others = 100-perc_unmodeled_fe-sum(value);
[~,b] = ismember(proteinid,gname_1);
proteinname = gname_2(b);

data = [others;value;perc_unmodeled_fe];

hold on;
b = bar(sI_res.k_cf(plotidx),transpose(data(:,plotidx)),'stacked');

for k = 1:length(b)
    b(k).FaceColor = color_set(k,:);
    b(k).FaceAlpha = 1;
    b(k).EdgeColor = 'w';
    b(k).EdgeAlpha = 0;
end

legend(['Others';proteinname;'Unmodeled'],'FontSize',6,'FontName','Helvetica','location','se');
legend('boxoff');
ylim([0 100]);
set(gca,'XColor','k');
set(gca,'YColor','k');
xlabel('Relative uptake','FontSize',6,'FontName','Helvetica');
ylabel('Iron usage fraction (%)','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
box off;



