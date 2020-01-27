%% addUnmodeledcofactor 
function model = addUnmodeledcofactor(model,cofactor_info,unmodeled_cofactor_id)

metlist = [cofactor_info.metid' {unmodeled_cofactor_id}];
coeflist = [(cofactor_info.mol_modeled-cofactor_info.mol_total)' 1];

rxnid = 'unmodeled_cofactor_formation';
model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);

