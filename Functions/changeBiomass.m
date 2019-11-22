%% addDummy 
function model = changeBiomass(model,f_modeled_protein,pseudo_biomass_rxn_id,protein_id)

metidx = ismember(model.mets,protein_id);
rxnidx = ismember(model.rxns,pseudo_biomass_rxn_id);
model.S(metidx,rxnidx) = full(model.S(metidx,rxnidx)) * (1-f_modeled_protein);
