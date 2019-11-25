%% updateComplexformation 
function model = updateComplexformation(model)

[~,~,raw] = xlsread('manual_update.xlsx','stoichiometry');

changed_enzymes = raw(2:end,1);
subunits = raw(2:end,2);
coeffs = raw(2:end,3);

for i = 1:length(changed_enzymes)
    enzymename = changed_enzymes{i};
    enzymeformationrxn = strcat(enzymename,'_formation');
    rxnidx = ismember(model.rxns,enzymeformationrxn);
    
    changed_subunits = subunits{i};
    changed_subunits = split(changed_subunits,'; ');
    changed_subunits = cellfun(@(x) strrep(x,'-','_'),changed_subunits,'UniformOutput',false);
    changed_subunits = cellfun(@(x) strcat(x,'_cofactorbound'),changed_subunits,'UniformOutput',false);
    
    changed_coeffs = coeffs{i};
    
    if length(changed_subunits) > 1
        changed_coeffs = split(changed_coeffs,'; ');
        for j = 1:length(changed_subunits)
            subunit_tmp = changed_subunits(j);
            coeff_tmp = -1 * str2double(changed_coeffs{j});
            metidx_tmp = ismember(model.mets,subunit_tmp);
            model.S(metidx_tmp,rxnidx) = coeff_tmp;
        end
    else
        coeff_tmp = -1 * changed_coeffs;
        metidx_tmp = ismember(model.mets,changed_subunits);
        model.S(metidx_tmp,rxnidx) = coeff_tmp;
    end
end
