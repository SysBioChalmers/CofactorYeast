load('sN1_fluxes.mat');
load('modelNoscapine.mat');
load('CofactorDataset.mat');
load('enzymedataNoscapine.mat');

[~,txt,~] = xlsread('Yeast8_Modification.xlsx','SGDgeneNames');
gname_1 = txt(2:end,1);
gname_2 = txt(2:end,2);
clear txt;
gname_2(ismember(gname_2,'')) = gname_1(ismember(gname_2,''));
gname_1 = [gname_1;{'Noscapine pathway'}];
gname_2 = [gname_2;{'Noscapine pathway'}];

flux = fluxes;
[~,b1] = size(flux);
mu = flux(strcmp(model.rxns,'r_2111'),:);
glc = -1*flux(strcmp(model.rxns,'r_1714'),:);
Noscapine = flux(strcmp(model.rxns,'new_r_eNoscapine'),:);
fe = -1*flux(ismember(model.rxnNames,'iron(3+) exchange'),:);

% extend CofactorDataset
[~,gprlst,~] = xlsread('New_pathway_noscapine.xlsx','Gpr_info');
[gpr_cfclst,~,~] = xlsread('New_pathway_noscapine.xlsx','Gpr_info');
gpr_gprlst = strtrim(gprlst(2:end,1));
gpr_cftlst = strtrim(gprlst(2:end,5));
gpr_cftsourcelst = strtrim(gprlst(2:end,7));
CofactorDataset.cofactor = [CofactorDataset.cofactor;gpr_cftlst(~ismember(gpr_cftlst,'-'))];
CofactorDataset.copy = [CofactorDataset.copy;gpr_cfclst(~ismember(gpr_cftlst,'-'))];
CofactorDataset.protein = [CofactorDataset.protein;gpr_gprlst(~ismember(gpr_cftlst,'-'))];
CofactorDataset.source = [CofactorDataset.source;gpr_cftsourcelst(~ismember(gpr_cftlst,'-'))];


conc_fe = zeros(length(model.genes),b1);
for j = 1:length(model.genes)
    protein_list_tmp = model.genes(j);
    conc_tmp = calculateCofactorUsage4protein(model,'FE',protein_list_tmp,CofactorDataset,flux);
    conc_fe(j,:) = conc_tmp;
end
tot_fe = fe./mu;
perc_fe = conc_fe./tot_fe*100;

new_idx = contains(model.genes,'uniprot_');
org_genes = model.genes(~new_idx);
org_genes_fe = perc_fe(~new_idx,:);
new_genes = {'Noscapine pathway'};
new_genes_fe = sum(perc_fe(new_idx,:));

data_genes = [org_genes;new_genes];
data_perc_fe = [org_genes_fe;new_genes_fe];

max_tmp = max(data_perc_fe(:,1),[],2);

top = 7;

[~,I] = sort(max_tmp,'descend');
proteinid = data_genes(I(1:top));
value = data_perc_fe(I(1:top),:);
others = 100-sum(value);
[~,b2] = ismember(proteinid,gname_1);
proteinname = gname_2(b2);

data = [others;value];

plotidx = 1:2:b1;
    
figure('Name','1');

color_set = [153,153,153
             228,26,28
             55,126,184
             152,78,163
             255,127,0
             77,175,74
             247,129,191
             166,86,40
             255,255,51]/255;

hold on;
b = bar(mu(plotidx),transpose(data(:,plotidx)),'stacked');

for k = 1:length(b)
    b(k).FaceColor = color_set(k,:);
    b(k).FaceAlpha = 0.8;
    b(k).EdgeColor = 'w';
    b(k).EdgeAlpha = 0;
end

legend(['Others';proteinname],'FontSize',6,'FontName','Helvetica','location','e');
legend('boxoff');
ylim([0 100]);
set(gca,'XColor','k');
set(gca,'YColor','k');
xlabel('Growth rate (/h)','FontSize',6,'FontName','Helvetica');
ylabel('Iron usage fraction (%)','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
box off;

set(gcf,'position',[0 0 440 140]);
set(gca,'position',[0.17 0.28 0.36 0.63]);


figure('Name','2');
subplot(3,1,1);
hold on;
box on;
plot(mu,Noscapine./glc,'-o','MarkerSize',2,'LineWidth',0.75,'Color',[55,126,184]/255);
xlim([0 0.4]);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
ylabel('Yield (mol/mol)','FontSize',7,'FontName','Helvetica');

subplot(3,1,2);
hold on;
box on;
plot(mu,Noscapine,'-o','MarkerSize',2,'LineWidth',0.75,'Color',[55,126,184]/255);
xlim([0 0.4]);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
ylabel(['Noscapine production flux',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');

subplot(3,1,3);
% load('sC_fluxes.mat');
% load('CofactorYeast.mat');
% mu_tmp = fluxes(strcmp(model.rxns,'r_2111'),:);
% fe_tmp = -1*fluxes(ismember(model.rxnNames,'iron(3+) exchange'),:);
hold on;
box on;
plot(mu,fe,'-o','MarkerSize',2,'LineWidth',0.75,'Color',[55,126,184]/255);
% plot(mu_tmp,fe_tmp,'-o','MarkerSize',2,'LineWidth',0.75,'Color',[135,135,135]/255);
xlim([0 0.4]);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
ylabel(['Iron uptake flux',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');

set(gcf,'position',[300 400 180 350]);