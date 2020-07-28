load('sB_res.mat');
load('CofactorYeast.mat');

protein_list = model.genes;

protein_conc = calculateProteinConc(model,protein_list,sB_res.fluxes);

ion_id_list = {'CA';'CU';'FE';'K';'MG';'MN';'NA';'ZN'};
ion_ex_list = {'r_4600'; ... % Ca(2+) exchange
               'r_4594'; ... % Cu(2+) exchange
               'r_1861'; ... % iron(2+) exchange
               'r_2020'; ... % potassium exchange
               'r_4597'; ... % Mg(2+) exchange
               'r_4595'; ... % Mn(2+) exchange
               'r_2049'; ... % sodium exchange
               'r_4596'};... % Zn(2+) exchange

mu = sB_res.fluxes(ismember(model.rxns,'r_2111'),:);

[~,b] = ismember(ion_ex_list,model.rxns);
ion_conc = -1*sB_res.fluxes(b,:)./mu;

cutoff_r = 0.9;

%% protein level

r_data = zeros(length(protein_list),length(ion_id_list));

for i = 1:length(ion_id_list)
    ion_conc_tmp = ion_conc(i,:);
    for j = 1:length(protein_list)
        r_tmp = corrcoef(ion_conc_tmp,protein_conc(j,:));
        r_data(j,i) = r_tmp(1,2);
    end
end

data_protein = protein_list(~all(isnan(r_data),2));
data_r_protein = r_data(~all(isnan(r_data),2),:);

[~,txt,~] = xlsread('Yeast8_Modification.xlsx','SGDgeneNames');
gname_1 = txt(2:end,1);
gname_2 = txt(2:end,2);
clear txt;
gname_2(ismember(gname_2,'')) = gname_1(ismember(gname_2,''));

[~,b] = ismember(data_protein,gname_1);
data_protein_id = gname_2(b);

data_protein_id_protein_level = data_protein_id(any(data_r_protein>cutoff_r,2));
data_r_protein_level = data_r_protein(any(data_r_protein>cutoff_r,2),:);

figure('Name','1');
h = heatmap(data_protein_id_protein_level,ion_id_list,data_r_protein_level','ColorMethod','count','CellLabelColor','none');
h.Colormap = hot;
h.Title = 'Pearson R';
set(h,'FontSize',8,'FontName','Helvetica');
set(gcf,'position',[500 300 1000 120]);
set(gca,'position',[0.02 0.4 0.93 0.5]);

%% pathway level
[~,txt1,~] = xlsread('kegg_pathway.xlsx','Sheet1');
[~,txt2,~] = xlsread('kegg_pathway.xlsx','Sheet2');
p2glist_p = txt1(:,1);
p2glist_g = txt1(:,2);
pname_p = txt2(:,1);
pname_name = txt2(:,2);
clear txt1 txt2;
global_pathway = {'path:sce01100';'path:sce01110';'path:sce01120';'path:sce01200';
'path:sce01210';'path:sce01212';'path:sce01230';'path:sce01220'};
idx_tmp = ~ismember(pname_p,global_pathway);
pname_p = pname_p(idx_tmp);
pname_name = pname_name(idx_tmp);
clear idx_tmp;

p2glist_g = cellfun(@(x) strrep(x,'sce:',''),p2glist_g,'UniformOutput',false);
pname_name = cellfun(@(x) strrep(x,' - Saccharomyces cerevisiae (budding yeast)',''),pname_name,'UniformOutput',false);

[uniprot_mw,txt,~] = xlsread('uniprot_protein_mw.xlsx');
uniprot_prot = txt(2:end,2);

[~,b] = size(sB_res.fluxes);
data_weight_pathway = zeros(length(pname_name),b); % g/gCDW
for i = 1:length(pname_name)
    sum_weight_tmp = zeros(1,b);
    genes_tmp = p2glist_g(ismember(p2glist_p,pname_p(i)));
    for j = 1:length(genes_tmp)
        if ismember(genes_tmp(j),protein_list)
            conc_tmp = protein_conc(ismember(protein_list,genes_tmp(j)),:);
            mw_tmp = uniprot_mw(ismember(uniprot_prot,genes_tmp(j)))/1000;
            sum_weight_tmp = sum_weight_tmp + conc_tmp*mw_tmp;
        end
    end
    data_weight_pathway(i,:) = sum_weight_tmp;
end

r_data = zeros(length(pname_name),length(ion_id_list));

for i = 1:length(ion_id_list)
    ion_conc_tmp = ion_conc(i,:);
    for j = 1:length(pname_name)
        r_tmp = corrcoef(ion_conc_tmp,data_weight_pathway(j,:));
        r_data(j,i) = r_tmp(1,2);
    end
end

data_pathway = pname_name(~all(isnan(r_data),2));
data_r_pathway = r_data(~all(isnan(r_data),2),:);

data_pathway = data_pathway(any(data_r_pathway>cutoff_r,2));
data_r_pathway = data_r_pathway(any(data_r_pathway>cutoff_r,2),:);

figure('Name','3');
h = heatmap(ion_id_list,data_pathway,data_r_pathway,'ColorMethod','count','CellLabelColor','none');
h.Colormap = hot;
h.Title = 'Pearson R';
set(h,'FontSize',6,'FontName','Helvetica');
set(gcf,'position',[200 200 500 700]);
set(gca,'position',[0.3 0.4 0.4 0.5]);











