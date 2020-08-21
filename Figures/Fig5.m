load('sN2_res.mat');
load('modelNoscapine.mat');

% selected_gene = {'Noscapine pathway';'Lys4';'Ole1';'Leu1';'Ilv3';'Cyc1';'Cob';'Rip1'};
% selected_geneID = {'Noscapine pathway';'YDR234W';'YGL055W';'YGL009C';'YJR016C';'YJR048W';'Q0105';'YEL024W'};
selected_gene = {'Noscapine pathway';'Lys4';'Ole1';'Leu1';'Ilv3';'ETC'};
selected_geneID = {'Noscapine pathway';'YDR234W';'YGL055W';'YGL009C';'YJR016C';'ETC'};

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
surf(-q_fe_list,data_mu,data_Noscapine,'EdgeColor','none');
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
set(gcf,'position',[0 300 150 130]);
set(gca,'position',[0.25 0.25 0.55 0.65]);

figure('Name','1-2d');
surf(-q_fe_list,data_mu,data_Noscapine,'EdgeColor','none');
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
set(gcf,'position',[300 300 150 100]);
set(gca,'position',[0.15 0.3 0.45 0.65]);

% fixed iron uptake rate
% fixed_fe = q_fe_list(floor(length(q_fe_list)/2));
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
set(gcf,'position',[600 300 150 130]);
set(gca,'position',[0.25 0.25 0.55 0.65]);

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
set(gcf,'position',[600 600 150 130]);
set(gca,'position',[0.25 0.25 0.55 0.65]);

%% fixed
[~,etc,~] = xlsread('kegg_pathway','ETC');
etc = unique(etc);

% for fixed fe uptake
fluxes_fixed_fe = sN2_res.fluxes(:,sN2_res.record == fixed_fe);
fe_fixed_fe = -fluxes_fixed_fe(ismember(model.rxnNames,'iron(3+) exchange'),:);
mu_fixed_fe = fluxes_fixed_fe(strcmp(model.rxns,'r_2111'),:);
Noscapine_fixed_fe = fluxes_fixed_fe(strcmp(model.rxns,'new_r_eNoscapine'),:);

% for fixed mu
fluxes_fixed_mu = sN2_res.fluxes(:,round(sN2_res.fluxes(strcmp(model.rxns,'r_2111'),:),3) == fixed_mu);
fe_fixed_mu = -fluxes_fixed_mu(ismember(model.rxnNames,'iron(3+) exchange'),:);
ub_fe_fixed_mu = -sN2_res.record(:,round(sN2_res.fluxes(strcmp(model.rxns,'r_2111'),:),3) == fixed_mu);
mu_fixed_mu = fluxes_fixed_mu(strcmp(model.rxns,'r_2111'),:);
Noscapine_fixed_mu = fluxes_fixed_mu(strcmp(model.rxns,'new_r_eNoscapine'),:);

figure('Name','2d_fixed_fe');
plot(mu_fixed_fe,Noscapine_fixed_fe,'-o','MarkerSize',0.5,'LineWidth',0.75,'Color','k');
xlim([0 0.4]);
ylim([0 0.1]);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate','FontSize',7,'FontName','Helvetica');
ylabel('Noscapine production','FontSize',7,'FontName','Helvetica');
title('Fixed iron UB','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 300 150 60]);
set(gca,'position',[0.17 0.28 0.3 0.63]);

figure('Name','2d_fixed_mu');
plot(ub_fe_fixed_mu,Noscapine_fixed_mu,'-o','MarkerSize',0.5,'LineWidth',0.75,'Color','k');
xlim([0 4e-4]);
ylim([0 0.1]);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Iron uptake UB','FontSize',7,'FontName','Helvetica');
ylabel('Noscapine production','FontSize',7,'FontName','Helvetica');
title('Fixed growth rate','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 0 150 60]);
set(gca,'position',[0.17 0.28 0.3 0.63]);


fluxes_combined = [fluxes_fixed_fe fluxes_fixed_mu];

conc_fe = calculateCofactorUsage4protein(model,'FE',model.genes,fluxes_combined);
tot_fe = [fe_fixed_fe fe_fixed_mu]./[mu_fixed_fe mu_fixed_mu];
perc_fe = conc_fe./tot_fe*100;
new_idx = contains(model.genes,'uniprot_');
org_genes = model.genes(~new_idx);
org_genes_fe = perc_fe(~new_idx,:);
new_genes = {'Noscapine pathway'};
new_genes_fe = sum(perc_fe(new_idx,:));

etc_idx = ismember(org_genes,etc);
org_genes_new = org_genes(~new_idx);
org_genes_fe_new = org_genes_fe(~new_idx,:);
etc_genes = {'ETC'};
etc_genes_fe = sum(org_genes_fe(etc_idx,:));

data_genes = [org_genes_new;new_genes;etc_genes];
data_perc_fe = [org_genes_fe_new;new_genes_fe;etc_genes_fe];




[~,b] = ismember(selected_geneID,data_genes);
value = data_perc_fe(b,:);
others = 100-sum(value);


color_set = [228,26,28
             55,126,184
             152,78,163
             255,127,0
             77,175,74
             247,129,191
             153,153,153]/255;


figure('Name','bar1');
hold on;
b = bar(mu_fixed_fe,transpose([value(:,1:length(mu_fixed_fe));others(:,1:length(mu_fixed_fe))]),'stacked');
for k = 1:length(b)
    b(k).FaceColor = color_set(k,:);
    b(k).FaceAlpha = 0.9;
    b(k).EdgeColor = 'w';
    b(k).EdgeAlpha = 0;
end
legend([selected_gene;'Others'],'FontSize',6,'FontName','Helvetica','location','e');
legend('boxoff');
ylim([0 100]);
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
ylabel('Iron usage fraction (%)','FontSize',7,'FontName','Helvetica');
title('Fixed UB of iron uptake','FontSize',7,'FontName','Helvetica');

box off;
set(gcf,'position',[600 300 440 140]);
set(gca,'position',[0.17 0.2 0.32 0.67]);


figure('Name','bar2');
hold on;
b = bar(ub_fe_fixed_mu,transpose([value(:,length(mu_fixed_fe)+1:end);others(:,length(mu_fixed_fe)+1:end)]),'stacked');
for k = 1:length(b)
    b(k).FaceColor = color_set(k,:);
    b(k).FaceAlpha = 0.9;
    b(k).EdgeColor = 'w';
    b(k).EdgeAlpha = 0;
end
legend([selected_gene;'Others'],'FontSize',6,'FontName','Helvetica','location','e');
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
set(gca,'position',[0.17 0.2 0.32 0.67]);




