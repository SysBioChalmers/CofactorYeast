load('sI_res.mat');
load('CofactorYeastExpand.mat');
load('enzymedataExpand.mat');

%% plot fluxes
fluxes = [sI_res.fluxes sI_res.flux_ref];
glc = -1*fluxes(strcmp(model.rxns,'r_1714'),:);
etoh = fluxes(strcmp(model.rxns,'r_1761'),:);
mu = fluxes(strcmp(model.rxns,'r_2111'),:);

figure('Name','1');
label_tmp = num2str(sI_res.k_cf);
label_tmp = strsplit(label_tmp);
label = [label_tmp 'Ref'];
subplot(3,1,1);
b1 = bar(glc,'FaceColor','flat','LineWidth',1);
b1.CData(end,:) = [1 1 1];
xticklabels(label);
ylim([15 20]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Glucose (mmol/gCDW/h)','FontSize',7,'FontName','Helvetica');
subplot(3,1,2);
b2 = bar(etoh,'FaceColor','flat','LineWidth',1);
b2.CData(end,:) = [1 1 1];
xticklabels(label);
ylim([23 33]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Ethanol (mmol/gCDW/h)','FontSize',7,'FontName','Helvetica');
subplot(3,1,3);
b3 = bar(mu,'FaceColor','flat','LineWidth',1);
b3.CData(end,:) = [1 1 1];
xticklabels(label);
ylim([0.39 0.42]);
yticks(0.39:0.01:0.42);
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

figure('Name','2');
data = log2(data_rel);
maxvalue = max(max(data));
minvalue = min(min(data));
data(data == -inf) = -maxvalue;

h = heatmap(label_tmp,protein_id_list,data,'Colormap',jet,'ColorMethod','count','CellLabelColor','none');
h.Title = 'log2FC';
set(h,'FontSize',6,'FontName','Helvetica');
% set(gcf,'position',[500 300 700 120]);
% set(gca,'position',[0.03 0.4 0.9 0.5]);

