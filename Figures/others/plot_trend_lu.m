load('sLU_res.mat');
load('CofactorYeast.mat');

ion = 'MG';

idx = contains(sLU_res.labels,ion);
fluxes = sLU_res.fluxes(:,idx);
labels = sLU_res.labels(1,idx);
labels = cellfun(@(x) x(strfind(x,'_')+1:end),labels,'UniformOutput',false);
labels = strrep(labels,'_','.');
labels = cellfun(@(x) str2double(x),labels,'UniformOutput',false);
lower_values = cell2mat(labels);

lower_values = lower_values(1,1:end-1);
fluxes = fluxes(:,1:end-1);
protein_conc = calculateProteinConc(model,model.genes,fluxes);

% remove low absolute protein level in reference
cutoff_low_abs = 0.05;
data_1 = protein_conc(:,1);
low_abs_value = quantile(data_1(data_1>0),cutoff_low_abs);
protein_list = model.genes(data_1>low_abs_value);
data_abs = protein_conc(data_1>low_abs_value,:);

% calculate relative value
data_rel = data_abs ./ data_abs(:,1);

% remove NAN data
data_rel_1 = data_rel(:,1);
protein_list = protein_list(~isnan(data_rel_1));
data_rel = data_rel(~isnan(data_rel_1),:);

% remove high relative value
cutoff_high_rel = 2;
highvalue_idx = any(data_rel>cutoff_high_rel,2);
protein_list = protein_list(~highvalue_idx);
data_rel = data_rel(~highvalue_idx,:);

[a,~] = size(data_rel);
for i = 1:a
    hold on;
    plot(lower_values,data_rel(i,:));
    
    
end

[~,txt,~] = xlsread('Yeast8_Modification.xlsx','SGDgeneNames');
gname_1 = txt(2:end,1);
gname_2 = txt(2:end,2);
clear txt;
gname_2(ismember(gname_2,'')) = gname_1(ismember(gname_2,''));

[~,b] = ismember(protein_list,gname_1);
protein_id_list = gname_2(b);



