%% calculateCompProteome 
function proteome_weight = calculateCompProteome(model,enzymedata,comps_id,flux)
% Unit: g/gCDW.

mu = flux(strcmp(model.rxns,'r_2111'),:);

[~, reactions] = collectCompartment(model,comps_id);

[~, b] = size(flux);

proteome_weight = zeros(1,b);
for i = 1:length(reactions)
    rxn = reactions{i};
    ezm = strcat(rxn,'_enzyme');
    dilrxn = strcat(ezm,'_dilution');
    conc = flux(ismember(model.rxns,dilrxn),:) ./ mu;
    mw = enzymedata.enzyme_MW(ismember(enzymedata.enzyme,ezm))/1000;
    proteome_weight = proteome_weight + mw * conc;
end








