%% calculateProteinConc 
function protein_conc = calculateProteinConc(model,protein_list,flux)
% Unit: mmol/gCDW.
% protein_list should be a column cell.
% The result shows the concentration of monomer.

[~, b] = size(flux);
mu = flux(ismember(model.rxns,'r_2111'),:);

protein_conc = zeros(length(protein_list),b);

for i = 1:length(protein_list)
    prot_tmp = protein_list(i);
    tsl_rxn_name = strcat('r_',strrep(prot_tmp,'-','_'),'_translated');
    tsl_flux_tmp = flux(ismember(model.rxns,tsl_rxn_name),:);
    protein_conc(i,:) = tsl_flux_tmp./mu;
end















