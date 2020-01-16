%% CofactorDataAnalysis
%  analysis of the collected cofactor dataset


% Timing: ~  s

tic;

%% Load data
load('CofactorDataset.mat');
load('Yeast8.mat');

%% Cofactor overview

overview = struct;
% Total types of cofactors detected in yeast
overview.cofactor = unique(CofactorDataset.cofactor);
% Frequency of each cofactor
overview.frequency = zeros(length(overview.cofactor),1);
for i = 1:length(overview.cofactor)
    tmp_cofactor = overview.cofactor(i);
    idx = ismember(CofactorDataset.cofactor,tmp_cofactor);
    protein_list = CofactorDataset.protein(idx);
    protein_list = unique(protein_list);
    overview.frequency(i) = length(protein_list);
end


idx_M = ismember(CofactorDataset.protein,org_model.genes);
M_cofactor = CofactorDataset.cofactor(idx_M);
M_protein = CofactorDataset.protein(idx_M);

overview_metabolism = struct;
% Total types of cofactors detected in yeast metabolism
overview_metabolism.cofactor = unique(M_cofactor);
% Frequency of each cofactor
overview_metabolism.frequency = zeros(length(overview_metabolism.cofactor),1);
for i = 1:length(overview_metabolism.cofactor)
    tmp_cofactor = overview_metabolism.cofactor(i);
    idx = ismember(M_cofactor,tmp_cofactor);
    protein_list = M_protein(idx);
    protein_list = unique(protein_list);
    overview_metabolism.frequency(i) = length(protein_list);
end




%% Cofactor abundance



[num,txt,~] = xlsread('Data.xlsx','Sheet6');
proteinlist = txt(2:end,1);
abundancelist = num;
clear num txt;

idx_mg = ismember(CofactorDataset.cofactor,'MG');
proteins = CofactorDataset.protein(idx_mg);
copies = CofactorDataset.copy(idx_mg);

abundance_mg = 0; % molecules per cell
for i = 1:length(proteins)
    proteinid = proteins(i);
    if ismember(proteinid,proteinlist)
        abund_tmp = abundancelist(ismember(proteinlist,proteinid));
        abundance_mg = abundance_mg + abund_tmp * copies(i);
    end
end
    

idx_M = ismember(CofactorDataset.protein,org_model.genes);
M_cofactor = CofactorDataset.cofactor(idx_M);
M_protein = CofactorDataset.protein(idx_M);
M_copy = CofactorDataset.copy(idx_M);

idx_mg = ismember(M_cofactor,'MG');
proteins = M_protein(idx_mg);
copies = M_copy(idx_mg);

abundance_mg = 0; % molecules per cell
for i = 1:length(proteins)
    proteinid = proteins(i);
    if ismember(proteinid,proteinlist)
        abund_tmp = abundancelist(ismember(proteinlist,proteinid));
        abundance_mg = abundance_mg + abund_tmp * copies(i);
    end
end






toc;


