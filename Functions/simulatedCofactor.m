%% simulatedCofactor 
function concentration = simulatedCofactor(model,cofactor_list,enzymedata,flux)
% Unit of concentration: mmol/gCDW.
% Cofactors are added to 1) each metabolic enzyme, 2) the biomass equation
% as unmodeled cofactor and 3) the dummy complex.

[~, n] = size(flux);

concentration = zeros(length(cofactor_list),n);
mu = flux(ismember(model.rxns,'r_2111'),:);

% 1) each metabolic enzyme
dil_rxns = model.rxns(contains(model.rxns,'_dilution'));

for i = 1:length(dil_rxns)
    enzymeid = strrep(dil_rxns(i),'_dilution','');
    flux_tmp = flux(ismember(model.rxns,dil_rxns(i)),:);
    
    idx_tmp = ismember(enzymedata.enzyme,enzymeid);
    type_tmp = enzymedata.cofactor_type(idx_tmp,:);
    copy_tmp = enzymedata.cofactor_copy(idx_tmp,:);
    type_tmp = type_tmp(copy_tmp ~= 0);
    copy_tmp = copy_tmp(copy_tmp ~= 0);
    conc_tmp = copy_tmp' * (flux_tmp ./ mu); %mmol/gCDW
    [~, b] = ismember(type_tmp,cofactor_list);
    concentration(b(b~=0),:) = concentration(b(b~=0),:) + conc_tmp(b~=0,:);
end

% 2) unmodeled cofactor
rxn_idx_tmp = ismember(model.rxns,'unmodeled_cofactor_formation');
met_idx_tmp = model.S(:,rxn_idx_tmp) < 0;
metlist = model.mets(met_idx_tmp);
slist = full(model.S(met_idx_tmp,rxn_idx_tmp));
flux_tmp = flux(ismember(model.rxns,'unmodeled_cofactor_formation'),:);

conc_tmp = -1 * slist * (flux_tmp ./ mu); %mmol/gCDW
[~, b] = ismember(metlist,cofactor_list);
concentration(b(b~=0),:) = concentration(b(b~=0),:) + conc_tmp(b~=0,:);

% 3) dummy cofactor

rxn_idx_tmp = ismember(model.rxns,'dilute_dummy_cofactor');
met_idx_tmp = model.S(:,rxn_idx_tmp) < 0;
metlist = model.mets(met_idx_tmp);
slist = full(model.S(met_idx_tmp,rxn_idx_tmp));
flux_tmp = flux(ismember(model.rxns,'dilute_dummy_cofactor'),:);

conc_tmp = -1 * slist * (flux_tmp ./ mu); %mmol/gCDW
[~, b] = ismember(metlist,cofactor_list);
concentration(b(b~=0),:) = concentration(b(b~=0),:) + conc_tmp(b~=0,:);






