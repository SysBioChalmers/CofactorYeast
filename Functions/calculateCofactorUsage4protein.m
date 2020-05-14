%% calculateCofactorUsage4protein 
function concentration = calculateCofactorUsage4protein(model,cofactor_type,protein_list,CofactorDataset,flux)
% Unit of concentration: mmol/gCDW.
% cofactor_type should be a string.
% protein_list should be a column cell.

if strcmp(cofactor_type,'CU')
    cofactor_type = {'CU_I';'CU_II'};
elseif strcmp(cofactor_type,'FE')
    cofactor_type = {'FE_III';'FE_II';'HEME_A';'HEME_C';'HEME_B';'ISC_2FE2S';'ISC_3FE4S';'ISC_4FE4S'};
else
    cofactor_type = {cofactor_type};
end

[~, n] = size(flux);

concentration = zeros(1,n);
mu = flux(ismember(model.rxns,'r_2111'),:);

for i = 1:length(protein_list)
    protein_id = protein_list{i};
    tansl_rxn_id = strcat('r_',strrep(protein_id,'-','_'),'_translated');
    if ismember(tansl_rxn_id,model.rxns)
        flux_tmp = flux(ismember(model.rxns,tansl_rxn_id),:);
        idx_tmp = ismember(CofactorDataset.protein,protein_id) & ismember(CofactorDataset.cofactor,cofactor_type);
        cofactor_tmp = CofactorDataset.cofactor(idx_tmp);
        copy_tmp = CofactorDataset.copy(idx_tmp);
        if ismember('ISC_2FE2S',cofactor_tmp)
            copy_tmp(ismember(cofactor_tmp,'ISC_2FE2S')) = copy_tmp(ismember(cofactor_tmp,'ISC_2FE2S'))*2;
        elseif ismember('ISC_3FE4S',cofactor_tmp)
            copy_tmp(ismember(cofactor_tmp,'ISC_3FE4S')) = copy_tmp(ismember(cofactor_tmp,'ISC_3FE4S'))*3;        
        elseif ismember('ISC_4FE4S',cofactor_tmp)
            copy_tmp(ismember(cofactor_tmp,'ISC_4FE4S')) = copy_tmp(ismember(cofactor_tmp,'ISC_4FE4S'))*4;        
        end
        conc_tmp = sum(copy_tmp) * (flux_tmp ./ mu); %mmol/gCDW
    else
        conc_tmp = 0;
    end
    concentration = concentration + conc_tmp;
end






