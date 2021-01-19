load('sLU_res.mat');
load('CofactorYeast.mat');

ion = 'FE';

idx = contains(sLU_res.labels,ion);
fluxes = sLU_res.fluxes(:,idx);

fluxes_ref = fluxes(:,1);

glc = -1*fluxes(strcmp(model.rxns,'r_1714'),:);
fluxes_low = fluxes(:,glc == max(glc));

protein_conc_ref = calculateProteinConc(model,model.genes,fluxes_ref);
protein_conc_low = calculateProteinConc(model,model.genes,fluxes_low);
protein_conc = [protein_conc_ref protein_conc_low];

% remove low absolute protein level in reference
cutoff_low_abs = 0.05;
low_abs_value = quantile(protein_conc_ref(protein_conc_ref>0),cutoff_low_abs);
protein_list = model.genes(protein_conc_ref>low_abs_value);
data_abs = protein_conc(protein_conc_ref>low_abs_value,:);

data_rel = data_abs(:,2)./data_abs(:,1);

[~,txt,~] = xlsread('Yeast8_Modification.xlsx','SGDgeneNames');
gname_1 = txt(2:end,1);
gname_2 = txt(2:end,2);
clear txt;
gname_2(ismember(gname_2,'')) = gname_1(ismember(gname_2,''));

[~,b] = ismember(protein_list,gname_1);
protein_id_list = gname_2(b);




