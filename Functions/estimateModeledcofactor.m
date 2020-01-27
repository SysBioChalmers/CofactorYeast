%% estimateModeledcofactor 
function cofactor_info = estimateModeledcofactor(model)
% Unit of abundance: molecule/cell.

[num, raw, ~] = xlsread('protein_abundance.xlsx','paxdb');
proteinlist = raw(2:end,2);
abundancelist = num(:,3); %unit:ppm
abundancelist = abundancelist*50*1e6/1e6; %unit:molecule/cell
% 50*1e6 protein molecules per cell in yeast (PMID: 24114984).

load('CofactorDataset.mat');

[~,raw,~] = xlsread('cofactorList.xlsx');
raw1 = raw(2:end,1);
raw2 = raw(2:end,2);
cofactor_list = raw1(~ismember(raw2,''));
cofactormetid = raw2(~ismember(raw2,''));

tot_abundance = zeros(length(cofactor_list),1);
modeled_abundance = zeros(length(cofactor_list),1);

for j = 1:length(cofactor_list)
    cofactorid = cofactor_list(j);
    
    idx_tmp = ismember(CofactorDataset.cofactor,cofactorid);
    proteins = CofactorDataset.protein(idx_tmp);
    copies = CofactorDataset.copy(idx_tmp);
    
    abundance_tot = 0; % molecule per cell
    abundance_modeled = 0; % molecule per cell
    for i = 1:length(proteins)
        proteinid = proteins(i);
        if ismember(proteinid,proteinlist)
            abund_tmp = abundancelist(ismember(proteinlist,proteinid));
            abundance_tot = abundance_tot + abund_tmp * copies(i);
            if ismember(proteinid,model.genes)
                abundance_modeled = abundance_modeled + abund_tmp * copies(i);
            end
        end
    end
    tot_abundance(j) = abundance_tot;
    modeled_abundance(j) = abundance_modeled;
end

% Convert molecule/cell to mmol/gCDW, 1 cell = 13 pg.
tot_mol = tot_abundance * 1e3 / (6.02e23*13e-12);
modeled_mol = modeled_abundance * 1e3 / (6.02e23*13e-12);

cofactor_info = struct();
cofactor_info.id = cofactor_list;
cofactor_info.metid = cofactormetid;
cofactor_info.abund_total = tot_abundance;
cofactor_info.abund_modeled = modeled_abundance;
cofactor_info.mol_total = tot_mol;
cofactor_info.mol_modeled = modeled_mol;