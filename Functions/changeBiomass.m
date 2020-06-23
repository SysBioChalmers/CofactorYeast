%% changeBiomass 
function model = changeBiomass(model,f_modeled_protein,pseudo_biomass_rxn_id,protein_id,unmodeled_cofactor_id,ion_id,cofactor_id)


rxnidx = ismember(model.rxns,pseudo_biomass_rxn_id);

% change protein coefficient
metidx = ismember(model.mets,protein_id);
model.S(metidx,rxnidx) = full(model.S(metidx,rxnidx)) * (1-f_modeled_protein);
% change ion coefficient
metidx = ismember(model.mets,ion_id);
model.S(metidx,rxnidx) = 0;
% add coefficient for unmodeled cofactor
metidx = ismember(model.mets,unmodeled_cofactor_id);
model.S(metidx,rxnidx) = -1;

% change cofactor coefficient
% metidx = ismember(model.mets,cofactor_id);
% model.S(metidx,rxnidx) = 0;
metidx = ismember(model.mets,cofactor_id);
model.S(ismember(model.metNames,'heme a [cytoplasm]'),model.S(metidx,:) > 0) = 0;




