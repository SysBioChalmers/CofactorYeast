%% Plot Biolog
load('sB_res.mat');

load('CofactorYeast.mat');
load('enzymedata.mat');

%% Compare exp data with sim data in terms of total metal ion abundances
% Sim data
ion_id_list = {'CA';'CU';'FE';'K';'MG';'MN';'NA';'ZN'};
ion_ex_list = {'r_4600'; ... % Ca(2+) exchange
               'r_4594'; ... % Cu(2+) exchange
               'r_1861'; ... % iron(2+) exchange
               'r_2020'; ... % potassium exchange
               'r_4597'; ... % Mg(2+) exchange
               'r_4595'; ... % Mn(2+) exchange
               'r_2049'; ... % sodium exchange
               'r_4596'};... % Zn(2+) exchange
sim_ion_abd = zeros(length(ion_id_list),length(sB_res.labels));
for i = 1:length(ion_id_list)
    idx_ion_ex_tmp = ismember(model.rxns,ion_ex_list(i));
    ion_uptake_tmp = -1*sB_res.fluxes(idx_ion_ex_tmp,:);
    tot_tmp = ion_uptake_tmp ./ sB_res.fluxes(strcmp(model.rxns,'r_2111'),:);
    % Convert mmol/gCDW to molecule/cell, 1 cell = 13 pg.
    sim_ion_abd(i,:) = tot_tmp * 6.02e23 * 13e-12 / 1e3;
end

% Exp dataset
exp_cofactor = cell(0,1);
exp_atomcell = zeros(0,1);
exp_source = cell(0,1);

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet1');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet2');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet3');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet4');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet5');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet6');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet7');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet8');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet9');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet10');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet11');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

clear atom_tmp idx_tmp cftrid_tmp tmp;

exp_ion_abd = nan(length(ion_id_list),length(unique(exp_source)));
for i = 1:length(ion_id_list)
    data_tmp = exp_atomcell(ismember(exp_cofactor,ion_id_list(i)))';
    exp_ion_abd(i,1:length(data_tmp)) = data_tmp;
end

% Plot
% figure('Name','1');
% hold on;
% h = boxplot(log10(sim_ion_abd)','Symbol','o','OutlierSize',3,'Widths',0.5,'Colors',[202,0,32]/255);
% set(h,{'linew'},{0.75});
% for i = 1:length(exp_cofactor)
%     x_tmp = find(ismember(ion_id_list,exp_cofactor{i}));
%     y_tmp = exp_atomcell(i);
%     scatter(x_tmp,log10(y_tmp),50,'x','LineWidth',1,'MarkerEdgeColor',[10,10,10]/255,'MarkerEdgeAlpha',0.8);
% end
% ylim([4.5 10.5]);
% xlim([0.1 length(ion_id_list)+0.9]);
% xticks(1:1:length(ion_id_list));
% set(gca, 'XColor','k');
% set(gca, 'YColor','k');
% set(gca,'XTickLabel',ion_id_list);
% set(gca,'FontSize',6,'FontName','Helvetica');
% ylabel('log10(atoms/cell)','FontSize',7,'FontName','Helvetica','Color','k');
% 
% set(gcf,'position',[500 100 220 150]);
% set(gca,'position',[0.11 0.1 0.87 0.85]);

figure('Name','1');
hold on;
h = boxplot(log10(sim_ion_abd)','Symbol','.','OutlierSize',10,'Widths',0.3,'Colors',[202,0,32]/255);
set(h,{'linew'},{0.5});
for i = 1:length(exp_cofactor)
    x_tmp = find(ismember(ion_id_list,exp_cofactor{i}))-0.3;
    y_tmp = exp_atomcell(i);
    scatter(x_tmp,log10(y_tmp),90,'.','LineWidth',1,'MarkerEdgeColor',[64,64,64]/255,'MarkerEdgeAlpha',0.8);
end
ylim([4.5 10.5]);
xlim([0.35 length(ion_id_list)+0.35]);
xticks(0.85:1:length(ion_id_list));
set(gca, 'XColor','k');
set(gca, 'YColor','k');
set(gca,'XTickLabel',ion_id_list);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('log10(atoms/cell)','FontSize',7,'FontName','Helvetica','Color','k');

set(gcf,'position',[500 100 220 130]);
set(gca,'position',[0.11 0.1 0.87 0.85]);