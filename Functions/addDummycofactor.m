%% addDummycofactor 
function model = addDummycofactor(model,cofactor_info)
% Assuming that the dummy complex has the same cofactor composition as
% metabolic protein pool.

rxnid = 'dilute_dummy_cofactor';

[~,idx] = ismember(cofactor_info.metid,model.metNames);

metlist = model.mets(idx)';
coeflist = -1*cofactor_info.mol_modeled';

model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);
