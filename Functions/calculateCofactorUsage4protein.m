%% calculateCofactorUsage4protein 
function concentration = calculateCofactorUsage4protein(model,cofactor_type,protein_list,flux)
% Unit of concentration: mmol/gCDW.
% cofactor_type should be a string.
% protein_list should be a column cell.

if strcmp(cofactor_type,'CU')
    cofactortype = {'CU_I';'CU_II'};
elseif strcmp(cofactor_type,'FE')
    cofactortype = {'FE_III';'FE_II';'HEME_A';'HEME_C';'HEME_B';'ISC_2FE2S';'ISC_3FE4S';'ISC_4FE4S'};
else
    cofactortype = {cofactor_type};
end

[~,raw,~] = xlsread('cofactorList.xlsx');
raw1 = raw(2:end,1);
raw2 = raw(2:end,2);
cofactorid = raw1(~ismember(raw2,''));
cofactormetname = raw2(~ismember(raw2,''));

cofactorname = cofactormetname(ismember(cofactorid,cofactortype));
cofactormetid = model.mets(ismember(model.metNames,cofactorname));

isc2fe = model.mets(ismember(model.metNames,cofactormetname(ismember(cofactorid,'ISC_2FE2S'))));
isc3fe = model.mets(ismember(model.metNames,cofactormetname(ismember(cofactorid,'ISC_3FE4S'))));
isc4fe = model.mets(ismember(model.metNames,cofactormetname(ismember(cofactorid,'ISC_4FE4S'))));

[~, n] = size(flux);

concentration = zeros(length(protein_list),n);
mu = flux(ismember(model.rxns,'r_2111'),:);

for i = 1:length(protein_list)
    
    disp(['Calculating usage for proteins ' cofactor_type ' : ' num2str(i) '/' num2str(length(protein_list))]);
    
    conc_tmp = zeros(1,n);
    protein_id = protein_list{i};
    
    tsl_met = strcat(strrep(protein_id,'-','_'),'_translated');
    
    cfbinding_rxn_idx = find(model.S(ismember(model.mets,tsl_met),:)<0);
    
    for j = 1:length(cfbinding_rxn_idx)
        flux_tmp = flux(cfbinding_rxn_idx(j),:);
        cf_met_idx = find(model.S(:,cfbinding_rxn_idx(j))<0);
        cf_list = model.mets(cf_met_idx);
        cf_copy_list = -1*full(model.S(cf_met_idx,cfbinding_rxn_idx(j)));
        
        if any(ismember(cf_list,cofactormetid))
            cf_list_new = cf_list(ismember(cf_list,cofactormetid));
            cf_copy_list_new = cf_copy_list(ismember(cf_list,cofactormetid));
            
            cf_copy_list_new(ismember(cf_list_new,isc2fe)) = cf_copy_list_new(ismember(cf_list_new,isc2fe)) * 2;
            cf_copy_list_new(ismember(cf_list_new,isc3fe)) = cf_copy_list_new(ismember(cf_list_new,isc3fe)) * 3;
            cf_copy_list_new(ismember(cf_list_new,isc4fe)) = cf_copy_list_new(ismember(cf_list_new,isc4fe)) * 4;
            
            conc_tmp = conc_tmp + sum(cf_copy_list_new) * (flux_tmp ./ mu); %mmol/gCDW
        end
        
    end
    concentration(i,:) = conc_tmp;
end
    
    
    
    
    
    
%     cofactorbinding_rxn_id = strcat('r_',strrep(protein_id,'-','_'),'_cofactorbound');
%     if ismember(cofactorbinding_rxn_id,model.rxns)
%         flux_tmp = flux(ismember(model.rxns,cofactorbinding_rxn_id),:);
%         idx_tmp = ismember(CofactorDataset.protein,protein_id) & ismember(CofactorDataset.cofactor,cofactor_type);
%         cofactor_tmp = CofactorDataset.cofactor(idx_tmp);
%         copy_tmp = CofactorDataset.copy(idx_tmp);
%         if ismember('ISC_2FE2S',cofactor_tmp)
%             copy_tmp(ismember(cofactor_tmp,'ISC_2FE2S')) = copy_tmp(ismember(cofactor_tmp,'ISC_2FE2S'))*2;
%         elseif ismember('ISC_3FE4S',cofactor_tmp)
%             copy_tmp(ismember(cofactor_tmp,'ISC_3FE4S')) = copy_tmp(ismember(cofactor_tmp,'ISC_3FE4S'))*3;        
%         elseif ismember('ISC_4FE4S',cofactor_tmp)
%             copy_tmp(ismember(cofactor_tmp,'ISC_4FE4S')) = copy_tmp(ismember(cofactor_tmp,'ISC_4FE4S'))*4;        
%         end
%         conc_tmp = conc_tmp + sum(copy_tmp) * (flux_tmp ./ mu); %mmol/gCDW
%     end
%     cofactorbindingwo_rxn_id = strcat('r_',strrep(protein_id,'-','_'),'_cofactorbound_withoutcofactor');
%     if ismember(cofactorbindingwo_rxn_id,model.rxns)
%         flux_tmp = flux(ismember(model.rxns,cofactorbindingwo_rxn_id),:);
%         idx_tmp = ismember(CofactorDataset.protein,protein_id) & ismember(CofactorDataset.cofactor,cofactor_type);
%         cofactor_tmp = CofactorDataset.cofactor(idx_tmp);
%         copy_tmp = CofactorDataset.copy(idx_tmp);
%         if ismember('ISC_2FE2S',cofactor_tmp)
%             copy_tmp(ismember(cofactor_tmp,'ISC_2FE2S')) = copy_tmp(ismember(cofactor_tmp,'ISC_2FE2S'))*2;
%         elseif ismember('ISC_3FE4S',cofactor_tmp)
%             copy_tmp(ismember(cofactor_tmp,'ISC_3FE4S')) = copy_tmp(ismember(cofactor_tmp,'ISC_3FE4S'))*3;        
%         elseif ismember('ISC_4FE4S',cofactor_tmp)
%             copy_tmp(ismember(cofactor_tmp,'ISC_4FE4S')) = copy_tmp(ismember(cofactor_tmp,'ISC_4FE4S'))*4;        
%         end
%         conc_tmp = conc_tmp + sum(copy_tmp) * (flux_tmp ./ mu); %mmol/gCDW
%     end
%     
%     concentration(i,:) = conc_tmp;







