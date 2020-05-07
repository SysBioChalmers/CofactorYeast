%% Plot Biolog
load('sB_res.mat');

load('CofactorYeast.mat');
load('enzymedata.mat');

%% plot
ion_id_list = {'K';'MG';'FE';'ZN';'CA';'MN';'CU';'NA'};
ion_ex_list = {'r_2020'; ... % potassium exchange
               'r_4597'; ... % Mg(2+) exchange
               'r_1861'; ... % iron(2+) exchange
               'r_4596'; ... % Zn(2+) exchange
               'r_4600'; ... % Ca(2+) exchange
               'r_4595'; ... % Mn(2+) exchange
               'r_4594'; ... % Cu(2+) exchange
               'r_2049'}; ... % sodium exchange

mu_list = sB_res.fluxes(ismember(model.rxns,'r_2111'),:);

source_label = cellfun(@(x) x(1:2),sB_res.labels,'UniformOutput',false);
C_idx = ismember(source_label,'C_');
N_idx = ismember(source_label,'N_');
P_idx = ismember(source_label,'P_');
S_idx = ismember(source_label,'S_');
C_clr = [228,26,28]/255;
N_clr = [55,126,184]/255;
P_clr = [255,127,0]/255;
S_clr = [152,78,163]/255;

% figure('Name','1');
% for i = 1:length(ion_ex_list)
%     ion_id = ion_id_list{i};
%     ion_exR = ion_ex_list{i};
%     ion_exR_flux = sB_res.fluxes(ismember(model.rxns,ion_exR),:);
%     ion_conc = -1 * ion_exR_flux ./ mu_list;
%     
%     ion_rel = log2(ion_conc ./ median(ion_conc));
%     
%     subplot(length(ion_id_list),1,i);
%     hold on;
%     box on;
%     scatter(1:1:length(sB_res.labels),ion_rel,10,mu_list,'filled');
%     h = colorbar;
%     colormap copper;
%     xlim([0 length(sB_res.labels)+1]);
%     ylim([-3 3]);
%     set(gca, 'XColor','k');
%     set(gca, 'YColor','k');
%     set(gca,'FontSize',6,'FontName','Helvetica');
% 	ylabel('Log2FC','FontSize',6,'FontName','Helvetica');
%     title(ion_id,'FontSize',7,'FontName','Helvetica','Color','k');
% end
% set(gcf,'position',[300 500 450 550]);

figure('Name','2');
for i = 1:length(ion_ex_list)
    ion_id = ion_id_list{i};
    ion_exR = ion_ex_list{i};
    ion_exR_flux = sB_res.fluxes(ismember(model.rxns,ion_exR),:);
    ion_conc = -1 * ion_exR_flux ./ mu_list;
    
    ion_rel = log2(ion_conc ./ median(ion_conc));
    
    subplot(length(ion_id_list)/2,2,i);
    hold on;
    box on;
    scatter(mu_list(C_idx),ion_rel(C_idx),10,'filled',...
        'MarkerEdgeColor','w','MarkerFaceColor',C_clr,...
        'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.6);
    scatter(mu_list(N_idx),ion_rel(N_idx),10,'filled',...
        'MarkerEdgeColor','w','MarkerFaceColor',N_clr,...
        'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.6);
    scatter(mu_list(P_idx),ion_rel(P_idx),10,'filled',...
        'MarkerEdgeColor','w','MarkerFaceColor',P_clr,...
        'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.6);
    scatter(mu_list(S_idx),ion_rel(S_idx),10,'filled',...
        'MarkerEdgeColor','w','MarkerFaceColor',S_clr,...
        'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.6);
    
    xlim([0 0.5]);
    ylim([-2.5 2.5]);
    set(gca, 'XColor','k');
    set(gca, 'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
    xlabel('Growth','FontSize',6,'FontName','Helvetica');
	ylabel('Log2FC','FontSize',6,'FontName','Helvetica');
    title(ion_id,'FontSize',7,'FontName','Helvetica','Color','k');
end
set(gcf,'position',[800 500 180 350]);