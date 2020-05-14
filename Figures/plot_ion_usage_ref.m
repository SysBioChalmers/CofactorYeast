%% plot ion distribution for reference

load('sR1.mat');
flux = fluxes;
load('CofactorYeast.mat');
load('CofactorDataset.mat')

[~,txt1,~] = xlsread('kegg_pathway.xlsx','Sheet1');
[~,txt2,~] = xlsread('kegg_pathway.xlsx','Sheet2');
p2glist_p = txt1(:,1);
p2glist_g = txt1(:,2);
pname_p = txt2(:,1);
pname_name = txt2(:,2);
clear txt1 txt2;
globel_pathway = {'path:sce01100';'path:sce01110';'path:sce01120';'path:sce01200';
'path:sce01210';'path:sce01212';'path:sce01230';'path:sce01220'};
idx_tmp = ~ismember(pname_p,globel_pathway);
pname_p = pname_p(idx_tmp);
pname_name = pname_name(idx_tmp);
clear idx_tmp;

p2glist_g = cellfun(@(x) strrep(x,'sce:',''),p2glist_g,'UniformOutput',false);
pname_name = cellfun(@(x) strrep(x,' - Saccharomyces cerevisiae (budding yeast)',''),pname_name,'UniformOutput',false);



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

top_proteins = 5;
top_pathways = 5;

color_proteins = [142,1,82]/255;
color_pathways = [39,100,25]/255;


for i = 1:length(ion_id_list)
    display([num2str(i),'/',num2str(length(ion_id_list))]);
    cofactor_type = ion_id_list{i};
    tot_tmp = -flux(ismember(model.rxns,ion_ex_list(i)))/flux(ismember(model.rxns,'r_2111'));
    
    % individual proteins
    conc_list_tmp = zeros(length(model.genes),1);
    for j = 1:length(model.genes)
        protein_list_tmp = model.genes(j);
        conc_tmp = calculateCofactorUsage4protein(model,cofactor_type,protein_list_tmp,CofactorDataset,flux);
        conc_list_tmp(j,1) = conc_tmp;
    end
    perc_list_tmp = conc_list_tmp/tot_tmp*100;
    [B,I] = sort(perc_list_tmp,'descend');
    value = B(1:top_proteins);
    proteinid = model.genes(I(1:top_proteins));
    [~,b] = ismember(proteinid,gname_1);
    proteinname = gname_2(b);
    proteinname(value == 0) = {'n/a'};
    
    subplot(2,length(ion_id_list),i);
    h1 = bar(1:length(proteinname),value,'FaceColor',color_proteins,'FaceAlpha',0.3,'EdgeColor',color_proteins,'LineWidth',0.5);
    set(gca,'XTick',1:1:length(proteinname));
    set(gca,'XTickLabel',proteinname);
    set(gca,'XColor','k');
    set(gca,'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
%     ylim([0.9 1.2]);
%     set(gca,'YTickLabel',0.9:0.1:1.2);
    ylabel('Fraction in total (%)','FontSize',6,'FontName','Helvetica','Color','k');
    title(cofactor_type,'FontSize',7,'FontName','Helvetica','Color','k');
    xtickangle(90);
    box off;
    
    % individual pathways
    conc_list_tmp = zeros(length(pname_p),1);
    for j = 1:length(pname_p)
        protein_list_tmp = p2glist_g(ismember(p2glist_p,pname_p(j)));
        conc_tmp = calculateCofactorUsage4protein(model,cofactor_type,protein_list_tmp,CofactorDataset,flux);
        conc_list_tmp(j,1) = conc_tmp;
    end
    perc_list_tmp = conc_list_tmp/tot_tmp*100;
    [B,I] = sort(perc_list_tmp,'descend');
    value = B(1:top_pathways);
    pathwayid = pname_name(I(1:top_pathways));
    pathwayid(value == 0) = {'n/a'};
    
    subplot(2,length(ion_id_list),i+length(ion_id_list));
    h2 = bar(1:length(pathwayid),value,'FaceColor',color_pathways,'FaceAlpha',0.3,'EdgeColor',color_pathways,'LineWidth',0.5);
    set(gca,'XTick',1:1:length(pathwayid));
    set(gca,'XTickLabel',pathwayid);
    set(gca,'XColor','k');
    set(gca,'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
%     ylim([0.9 1.2]);
%     set(gca,'YTickLabel',0.9:0.1:1.2);
    ylabel('Fraction in total (%)','FontSize',6,'FontName','Helvetica','Color','k');
    title(cofactor_type,'FontSize',7,'FontName','Helvetica','Color','k');
    xtickangle(90);
    box off;
end

set(gcf,'position',[300 500 650 380]);



