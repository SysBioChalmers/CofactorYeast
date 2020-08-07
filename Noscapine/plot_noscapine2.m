load('sN2_res.mat');
load('modelNoscapine.mat');

q_fe_list = unique(sN2_res.record);
tot_raw = sum(sN2_res.record == min(q_fe_list));
data_mu = zeros(tot_raw,length(q_fe_list));
data_Noscapine = zeros(tot_raw,length(q_fe_list));

for i = 1:length(q_fe_list)
    idx_tmp = sN2_res.record == q_fe_list(i);
    flux_tmp = sN2_res.fluxes(:,idx_tmp);
    mu = flux_tmp(strcmp(model.rxns,'r_2111'),:);
    glc = -1*flux_tmp(strcmp(model.rxns,'r_1714'),:);
    Noscapine = flux_tmp(strcmp(model.rxns,'new_r_eNoscapine'),:);
%     Noscapine = Noscapine./glc;
    data_mu(:,i) = [mu nan(1,tot_raw-length(mu))];
    data_Noscapine(:,i) = [Noscapine nan(1,tot_raw-length(Noscapine))];
end

%% 3d plot
figure('Name','1');
surf(-q_fe_list,data_mu,data_Noscapine);
hold on;
grid on;
az = -145;
el = 25;
view(az, el);
xlim([0 4e-4]);
zlim([0 0.4]);
zlim([0 0.1]);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel(['UB of iron uptake',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');
ylabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
zlabel(['Noscapine production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[0 600 250 200]);
set(gca,'position',[0.25 0.25 0.65 0.65]);

figure('Name','1-2d');
surf(-q_fe_list,data_mu,data_Noscapine);
hold on;
grid on;
colorbar;
az = -180;
el = 90;
view(az, el);
xlim([0 4e-4]);
zlim([0 0.4]);
zlim([0 0.1]);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel(['UB of iron uptake',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');
ylabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
zlabel(['Noscapine production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 600 250 200]);
set(gca,'position',[0.15 0.25 0.55 0.65]);

% fixed iron uptake rate
fixed_fe = q_fe_list(10);
figure('Name','2');
surf(-q_fe_list,data_mu,data_Noscapine,'FaceAlpha',0.6,'EdgeColor','none');
hold on;
grid on;
az = -145;
el = 25;
view(az, el);
xlim([0 4e-4]);
ylim([0 0.4]);
zlim([0 0.1]);
plot3(-fixed_fe*[1,1],[0,0],[0,0.1],'k-','LineWidth',1);
plot3(-fixed_fe*[1,1],[0.4,0.4],[0,0.1],'k-','LineWidth',1);
plot3(-fixed_fe*[1,1],[0.4,0],[0,0],'k-','LineWidth',1);
plot3(-fixed_fe*[1,1],[0.4,0],[0.1,0.1],'k-','LineWidth',1);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel(['UB of iron uptake',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');
ylabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
zlabel(['Noscapine production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[0 300 250 200]);
set(gca,'position',[0.25 0.25 0.65 0.65]);

% fixed growth rate
fixed_mu = 0.2;
figure('Name','3');
surf(-q_fe_list,data_mu,data_Noscapine,'FaceAlpha',0.6,'EdgeColor','none');
hold on;
grid on;
az = -145;
el = 25;
view(az, el);
xlim([0 4e-4]);
ylim([0 0.4]);
zlim([0 0.1]);
plot3([0,0],fixed_mu*[1,1],[0,0.1],'k-','LineWidth',1);
plot3([0,4e-4],fixed_mu*[1,1],[0,0],'k-','LineWidth',1);
plot3([4e-4,4e-4],fixed_mu*[1,1],[0,0.1],'k-','LineWidth',1);
plot3([4e-4,0],fixed_mu*[1,1],[0.1,0.1],'k-','LineWidth',1);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel(['UB of iron uptake',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');
ylabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
zlabel(['Noscapine production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[0 0 250 200]);
set(gca,'position',[0.25 0.25 0.65 0.65]);

%% 2d plot
load('CofactorDataset.mat');
load('enzymedataNoscapine.mat');

load('cofactor_info.mat');
idx_ion = ismember(cofactor_info.element_id,'FE');
tot_unmodeled_tmp = cofactor_info.element_abund_total(idx_ion) - cofactor_info.element_abund_modeled(idx_ion);
tot_unmodeled_tmp = tot_unmodeled_tmp * 1e3 / (6.02e23*13e-12);

[~,txt,~] = xlsread('Yeast8_Modification.xlsx','SGDgeneNames');
gname_1 = txt(2:end,1);
gname_2 = txt(2:end,2);
clear txt;
gname_2(ismember(gname_2,'')) = gname_1(ismember(gname_2,''));
gname_1 = [gname_1;{'Noscapine pathway'}];
gname_2 = [gname_2;{'Noscapine pathway'}];
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

color_set = [0,0,0
             153,153,153
             228,26,28
             55,126,184
             152,78,163
             255,127,0
             77,175,74
             247,129,191
             166,86,40
             255,255,51]/255;

% for fixed fe uptake
fluxes_fixed_fe = sN2_res.fluxes(:,sN2_res.record == fixed_fe);
fe_fixed_fe = -fluxes_fixed_fe(ismember(model.rxnNames,'iron(3+) exchange'),:);
mu_fixed_fe = fluxes_fixed_fe(strcmp(model.rxns,'r_2111'),:);
Noscapine_fixed_fe = fluxes_fixed_fe(strcmp(model.rxns,'new_r_eNoscapine'),:);

% 2dplot
figure('Name','2d_fixed_fe');
plot(mu_fixed_fe,Noscapine_fixed_fe,'-o','MarkerSize',2,'LineWidth',0.75,'Color','k');
xlim([0 0.4]);
ylim([0 0.1]);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
ylabel(['Noscapine production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');
title('Fixed iron UB','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[300 300 250 120]);
set(gca,'position',[0.17 0.28 0.36 0.63]);

conc_fe = calculateCofactorUsage4protein(model,'FE',model.genes,fluxes_fixed_fe);
tot_fe = fe_fixed_fe./mu_fixed_fe;
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
perc_unmodeled_fe = tot_unmodeled_tmp./tot_fe*100;
others = 100-perc_unmodeled_fe-sum(value);
[~,b2] = ismember(proteinid,gname_1);
proteinname = gname_2(b2);    
figure('Name','2d1');
hold on;
b = bar(mu_fixed_fe,transpose([perc_unmodeled_fe;others;value]),'stacked');
for k = 1:length(b)
    b(k).FaceColor = color_set(k,:);
    b(k).FaceAlpha = 0.8;
    b(k).EdgeColor = 'w';
    b(k).EdgeAlpha = 0;
end
legend(['Unmodeled';'Others';proteinname],'FontSize',6,'FontName','Helvetica','location','e');
legend('boxoff');
ylim([0 100]);
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
ylabel('Iron usage fraction (%)','FontSize',7,'FontName','Helvetica');
title('Fixed iron UB','FontSize',7,'FontName','Helvetica');


box off;
set(gcf,'position',[600 300 440 140]);
set(gca,'position',[0.17 0.28 0.36 0.63]);


% for fixed mu
fluxes_fixed_mu = sN2_res.fluxes(:,round(sN2_res.fluxes(strcmp(model.rxns,'r_2111'),:),3) == fixed_mu);
fe_fixed_mu = -fluxes_fixed_mu(ismember(model.rxnNames,'iron(3+) exchange'),:);
ub_fe_fixed_mu = -sN2_res.record(:,round(sN2_res.fluxes(strcmp(model.rxns,'r_2111'),:),3) == fixed_mu);
mu_fixed_mu = fluxes_fixed_mu(strcmp(model.rxns,'r_2111'),:);
Noscapine_fixed_mu = fluxes_fixed_mu(strcmp(model.rxns,'new_r_eNoscapine'),:);

% 2dplot
figure('Name','2d_fixed_mu');
plot(ub_fe_fixed_mu,Noscapine_fixed_mu,'-o','MarkerSize',2,'LineWidth',0.75,'Color','k');
xlim([0 4e-4]);
ylim([0 0.1]);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel(['UB of iron uptake',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');
ylabel(['Noscapine production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');
title('Fixed growth rate','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[300 0 250 120]);
set(gca,'position',[0.17 0.28 0.36 0.63]);


conc_fe = calculateCofactorUsage4protein(model,'FE',model.genes,fluxes_fixed_mu);
tot_fe = fe_fixed_mu./mu_fixed_mu;
perc_fe = conc_fe./tot_fe*100;
new_idx = contains(model.genes,'uniprot_');
org_genes = model.genes(~new_idx);
org_genes_fe = perc_fe(~new_idx,:);
new_genes = {'Noscapine pathway'};
new_genes_fe = sum(perc_fe(new_idx,:));
data_genes = [org_genes;new_genes];
data_perc_fe = [org_genes_fe;new_genes_fe];
max_tmp = max(data_perc_fe(:,end),[],2);
top = 7;
[~,I] = sort(max_tmp,'descend');
proteinid = data_genes(I(1:top));
value = data_perc_fe(I(1:top),:);
perc_unmodeled_fe = tot_unmodeled_tmp./tot_fe*100;
others = 100-perc_unmodeled_fe-sum(value);
[~,b2] = ismember(proteinid,gname_1);
proteinname = gname_2(b2);    
figure('Name','2d2');
hold on;
b = bar(ub_fe_fixed_mu,transpose([perc_unmodeled_fe;others;value]),'stacked');
for k = 1:length(b)
    b(k).FaceColor = color_set(k,:);
    b(k).FaceAlpha = 0.8;
    b(k).EdgeColor = 'w';
    b(k).EdgeAlpha = 0;
end
legend(['Unmodeled';'Others';proteinname],'FontSize',6,'FontName','Helvetica','location','e');
legend('boxoff');
ylim([0 100]);
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel(['UB of iron uptake',char(13,10)','(mmol/gCDW/h)'],'FontSize',7,'FontName','Helvetica');
ylabel('Iron usage fraction (%)','FontSize',7,'FontName','Helvetica');
title('Fixed growth rate','FontSize',7,'FontName','Helvetica');

box off;
set(gcf,'position',[600 0 440 140]);
set(gca,'position',[0.17 0.28 0.36 0.63]);




