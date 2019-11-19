%% Process protein stoichiometry data

% Determine protein stoichiometry for each protein based on the collected
% data from PDB database. Almost all the proteins have multiple PDB IDs and
% therefore multiple stoichiometric coefficients.
% If all the coefficients are the same value, then it would be selected as
% the final value for the protein; if a protein has various coefficients,
% then the mimimum would be selected.
% Besides, the proteins having various coefficients are collected in the
% file "Undetermined_stoichiometry.mat".

[num,txt,~] = xlsread('ProteinStoichiometry.xlsx','All_GEM_proteins');
protein_list = txt(2:end,1);
pdb_list = txt(2:end,2);
s_coefficent_list = num(1:end,3);
clear num txt;

unique_protein_list = unique(protein_list);

Determined_stoichiometry = struct;
Determined_stoichiometry.protein = cell(0,1);
Determined_stoichiometry.coefficient = zeros(0,1);

Undetermined_stoichiometry = struct;
Undetermined_stoichiometry.protein = cell(0,1);
Undetermined_stoichiometry.pdbid = cell(0,1);
Undetermined_stoichiometry.coefficient = zeros(0,1);

for i = 1:length(unique_protein_list)
    proteinname = unique_protein_list(i);
    idx = ismember(protein_list,proteinname);
    s_coefficents = s_coefficent_list(idx);
    proteins = protein_list(idx);
    pdbs = pdb_list(idx);
    if length(unique(s_coefficents)) == 1
        Determined_stoichiometry.protein = [Determined_stoichiometry.protein;proteinname];
        Determined_stoichiometry.coefficient = [Determined_stoichiometry.coefficient;unique(s_coefficents)];
    else
        Determined_stoichiometry.protein = [Determined_stoichiometry.protein;proteinname];
        Determined_stoichiometry.coefficient = [Determined_stoichiometry.coefficient;min(s_coefficents)];
        Undetermined_stoichiometry.protein = [Undetermined_stoichiometry.protein;proteins];
        Undetermined_stoichiometry.pdbid = [Undetermined_stoichiometry.pdbid;pdbs];
        Undetermined_stoichiometry.coefficient = [Undetermined_stoichiometry.coefficient;s_coefficents];
    end
end

save('Determined_stoichiometry.mat','Determined_stoichiometry');
save('Undetermined_stoichiometry.mat','Undetermined_stoichiometry');
clear;