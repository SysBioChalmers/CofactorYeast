%% plot ion distribution for lower uptake

load('sLU_res.mat');
load('CofactorYeast.mat');
load('CofactorDataset.mat')
load('enzymedata.mat');
load('cofactor_info.mat');

[~,txt,~] = xlsread('Yeast8_Modification.xlsx','SGDgeneNames');
gname_1 = txt(2:end,1);
gname_2 = txt(2:end,2);
clear txt;
gname_2(ismember(gname_2,'')) = gname_1(ismember(gname_2,''));


ion_id_list = {'CA';'CU';'FE';'K';'MG';'MN';'NA';'ZN'};
ion_ex_list = {'r_4600'; ... % Ca(2+) exchange
               'r_4594'; ... % Cu(2+) exchange
               'r_1861'; ... % iron(2+) exchange
               'r_2020'; ... % potassium exchange
               'r_4597'; ... % Mg(2+) exchange
               'r_4595'; ... % Mn(2+) exchange
               'r_2049'; ... % sodium exchange
               'r_4596'};... % Zn(2+) exchange


lbl = sLU_res.labels;
lbl = cellfun(@(x) x(1:strfind(x,'_')-1),lbl,'UniformOutput',false);
lbl = unique(lbl);

top_proteins = 6;

color_set = [153,153,153
             228,26,28
             55,126,184
             152,78,163
             255,127,0
             77,175,74
             247,129,191
             255,255,51
             166,86,40]/255;


for i = 1:length(lbl)
    ion = lbl{i};
    idx = contains(sLU_res.labels,ion);
    fluxes_tmp = sLU_res.fluxes(:,idx);
    labels_tmp = sLU_res.labels(1,idx);
    labels_tmp = cellfun(@(x) x(strfind(x,'_')+1:end),labels_tmp,'UniformOutput',false);
    labels_tmp = strrep(labels_tmp,'_','.');
    labels_tmp = cellfun(@(x) str2double(x),labels_tmp,'UniformOutput',false);
    lower_values = cell2mat(labels_tmp);
    
    display([num2str(i),'/',num2str(length(ion_id_list))]);
    
    tot_tmp = -fluxes_tmp(ismember(model.rxns,ion_ex_list(ismember(ion_id_list,ion))),:)./fluxes_tmp(ismember(model.rxns,'r_2111'),:);
    idx_ion = ismember(cofactor_info.element_id,ion);
    tot_unmodeled_tmp = cofactor_info.element_abund_total(idx_ion) - cofactor_info.element_abund_modeled(idx_ion);
    tot_unmodeled_tmp = tot_unmodeled_tmp * 1e3 / (6.02e23*13e-12);
    tot_modeled_tmp = tot_tmp - tot_unmodeled_tmp;
    
    % individual proteins
    conc_list_tmp = zeros(length(model.genes),length(lower_values));
    for j = 1:length(model.genes)
        protein_list_tmp = model.genes(j);
        conc_tmp = calculateCofactorUsage4protein(model,ion,protein_list_tmp,CofactorDataset,fluxes_tmp);
        conc_list_tmp(j,:) = conc_tmp;
    end
    perc_list_tmp = conc_list_tmp./tot_modeled_tmp*100;
    
    max_tmp = max(perc_list_tmp,[],2);
    
    if length(find(max_tmp)) >= top_proteins
        top = top_proteins;
    else
        top = length(find(perc_list_tmp(:,1)));
    end
    
    [~,I] = sort(max_tmp,'descend');
    proteinid = model.genes(I(1:top));
    value = perc_list_tmp(I(1:top),:);
    others = 100-sum(value);
    [~,b] = ismember(proteinid,gname_1);
    proteinname = gname_2(b);
    
    data = [others;value];
    
    plotidx = 1:2:length(lower_values);
    
    subplot(2,length(lbl),i);
    hold on;
    b = bar(lower_values(plotidx),transpose(data(:,plotidx)),'stacked');
    
    for k = 1:length(b)
        b(k).FaceColor = color_set(k,:);
        b(k).FaceAlpha = 0.8;
        b(k).EdgeColor = 'w';
        b(k).EdgeAlpha = 0;
    end
    
    legend(['Others';proteinname],'FontSize',6,'FontName','Helvetica','location','se');
    legend('boxoff');
    ylim([0 100]);
    set(gca,'XColor','k');
    set(gca,'YColor','k');
    xlabel('Relative uptake','FontSize',6,'FontName','Helvetica');
    if i == 1
        ylabel('Ion usage fraction (%)','FontSize',6,'FontName','Helvetica');
    end
    set(gca,'FontSize',6,'FontName','Helvetica');
    title(ion,'FontSize',7,'FontName','Helvetica','Color','k');
    box off;
end

set(gcf,'position',[300 500 650 180]);



