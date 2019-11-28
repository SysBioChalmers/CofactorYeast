%% calculateMW 
function tot_protein_cost = calculateProteinCost(model,fluxes,enzymedata)

non_zero_rxns = model.rxns(fluxes ~= 0);
non_zero_fluxes = fluxes(fluxes ~= 0);

tot_protein_cost = 0;

for i = 1:length(non_zero_rxns)
    rxnid = non_zero_rxns{i};
    flux = non_zero_fluxes(i);
    enzymeid = strcat(rxnid,'_enzyme');
    
    if ismember(enzymeid,enzymedata.enzyme)
        idx_tmp = ismember(enzymedata.enzyme,enzymeid);
        proteincost = enzymedata.enzyme_MW(idx_tmp) / enzymedata.kcat(idx_tmp);
        tot_protein_cost = tot_protein_cost + proteincost * flux;
    end
end

tot_protein_cost = tot_protein_cost / 1000; %unit:gProtein/gCDW
