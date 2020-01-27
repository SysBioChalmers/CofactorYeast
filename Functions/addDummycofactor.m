%% addDummycofactor 
function model = addDummycofactor(model,cofactor_info)
% Assuming that the dummy complex has the same cofactor composition as
% metabolic protein pool.

rxnid = 'dilute_dummy_cofactor';

metlist = cofactor_info.metid';
coeflist = -1*cofactor_info.mol_modeled';

model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);
