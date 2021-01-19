%% Plot Biolog
load('sB_res.mat');

load('CofactorYeast.mat');
load('enzymedata.mat');

%% plot
ion_id_list = {'CA';'CU';'FE';'K';'MG';'MN';'NA';'ZN'};
ion_ex_list = {'r_4600'; ... % Ca(2+) exchange
               'r_4594'; ... % Cu(2+) exchange
               'r_1861'; ... % iron(2+) exchange
               'r_2020'; ... % potassium exchange
               'r_4597'; ... % Mg(2+) exchange
               'r_4595'; ... % Mn(2+) exchange
               'r_2049'; ... % sodium exchange
               'r_4596'};... % Zn(2+) exchange



% reference glucose
glc_id = model.mets(ismember(model.metNames,'D-glucose [extracellular]'));
glc_id = strrep(glc_id,'[','_');
glc_id = strrep(glc_id,']','');
glc_id = strcat('C_',glc_id);
mu_glc = sB_res.fluxes(ismember(model.rxns,'r_2111'),ismember(sB_res.labels,glc_id));
ionEX_flux_glc = sB_res.fluxes(ismember(model.rxns,ion_ex_list),ismember(sB_res.labels,glc_id));
ion_conc_glc = -ionEX_flux_glc / mu_glc;

% re-organize dataset
mu_all = sB_res.fluxes(ismember(model.rxns,'r_2111'),:);
ionEX_flux_all = sB_res.fluxes(ismember(model.rxns,ion_ex_list),:);
ion_conc_all = -ionEX_flux_all ./ mu_all;
% order by source label
[label1,I1] = sort(sB_res.labels);
mu_all1 = mu_all(:,I1);
ion_conc_all1 = ion_conc_all(:,I1);
% order by growth rate
[mu_all2,I2] = sort(mu_all1);
label2 = label1(:,I2);
ion_conc_all2 = ion_conc_all1(:,I2);

mu_list = mu_all2;
label_list = label2;
ion_conc_list = ion_conc_all2;
ref_ion_conc_list = ion_conc_list ./ ion_conc_glc;
ref_mu_list = mu_list / mu_glc;
log2_ref_ion_conc_list = log2(ref_ion_conc_list);
log2_ref_mu_list = log2(ref_mu_list);


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

source_label = cellfun(@(x) x(1:2),sB_res.labels,'UniformOutput',false);
C_idx = ismember(source_label,'C_');
N_idx = ismember(source_label,'N_');
P_idx = ismember(source_label,'P_');
S_idx = ismember(source_label,'S_');
C_clr = [228,26,28]/255;
N_clr = [55,126,184]/255;
P_clr = [255,127,0]/255;
S_clr = [152,78,163]/255;

for i = 1:length(ion_id_list)
    ion_rel = log2_ref_ion_conc_list(i,:);
    
    subplot(1,length(ion_id_list),i);
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
    ylim([-2.2 2.2]);
    set(gca, 'XColor','k');
    set(gca, 'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
    xlabel('Growth (/h)','FontSize',6,'FontName','Helvetica');
    if i == 1
        ylabel('Log2FC','FontSize',6,'FontName','Helvetica');
    end
    title(ion_id_list{i},'FontSize',7,'FontName','Helvetica','Color','k');
end
set(gcf,'position',[800 500 550 80]);

% heatmap
% colormap
data = [log2_ref_ion_conc_list;log2_ref_mu_list];

maxclr = [103,0,31]/255;
mdlclr = [255,255,255]/255;
minclr = [5,48,97]/255;
tmp11 = linspace(minclr(1),mdlclr(1),128)';
tmp12 = linspace(mdlclr(1),maxclr(1),128)';
tmp21 = linspace(minclr(2),mdlclr(2),128)';
tmp22 = linspace(mdlclr(2),maxclr(2),128)';
tmp31 = linspace(minclr(3),mdlclr(3),128)';
tmp32 = linspace(mdlclr(3),maxclr(3),128)';
tmp1 = [tmp11;tmp12(2:end)];
tmp2 = [tmp21;tmp22(2:end)];
tmp3 = [tmp31;tmp32(2:end)];
clrmap = [tmp1 tmp2 tmp3];

maxvalue = max(max(data));
minvalue = min(min(data));
minposition = 1;
mdlposition = find(any((clrmap == mdlclr)'));
maxposition = length(clrmap);
if abs(minvalue) >= abs(maxvalue)
    startposition = minposition;
    endposition = round(abs(maxvalue/minvalue)*(maxposition-mdlposition))+mdlposition;
else
    endposition = maxposition;
    startposition = mdlposition-floor(abs(minvalue/maxvalue)*mdlposition);
end
clrmap = clrmap(startposition:endposition,:);

figure('Name','3');
h = heatmap(label_list,ion_id_list,log2_ref_ion_conc_list,'Colormap',clrmap,'ColorMethod','count','CellLabelColor','none');
h.Title = 'log2FC';
set(h,'FontSize',6,'FontName','Helvetica');
set(gcf,'position',[500 300 700 120]);
set(gca,'position',[0.03 0.4 0.9 0.5]);

