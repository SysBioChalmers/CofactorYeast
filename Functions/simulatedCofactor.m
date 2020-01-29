%% simulatedCofactor 
function [list, concentration] = simulatedCofactor(model,enzymedata,flux)
% Unit of concentration: mmol/gCDW.
% Cofactors are added to 1) each metabolic enzyme, 2) the biomass equation
% as unmodeled cofactor and 3) the dummy complex.

list = cell(0,1);
concentration = zeros(0,1);

mu = flux(ismember(model.rxns,'r_2111'));

% 1) each metabolic enzyme
dil_rxns = model.rxns(contains(model.rxns,'_dilution'));

for i = 1:length(dil_rxns)
    enzymeid = strrep(dil_rxns(i),'_dilution','');
    flux_tmp = flux(ismember(model.rxns,dil_rxns(i)));
    if flux_tmp ~= 0
        idx_tmp = ismember(enzymedata.enzyme,enzymeid);
        type_tmp = enzymedata.cofactor_type(idx_tmp,:);
        copy_tmp = enzymedata.cofactor_copy(idx_tmp,:);
        type_tmp = type_tmp(copy_tmp ~= 0);
        copy_tmp = copy_tmp(copy_tmp ~= 0);
        if ~isempty(type_tmp)
            conc_tmp = copy_tmp * flux_tmp / mu; %mmol/gCDW
            for j = 1:length(type_tmp)
                if ~ismember(type_tmp(j),list)
                    list = [list;type_tmp(j)];
                    concentration = [concentration;conc_tmp(j)];
                else
                    idx_tmptmp = ismember(list,type_tmp(j));
                    concentration(idx_tmptmp) = concentration(idx_tmptmp) + conc_tmp(j);
                end
            end
        end
    end
end

% 2) unmodeled cofactor

rxn_idx_tmp = ismember(model.rxns,'unmodeled_cofactor_formation');
met_idx_tmp = model.S(:,rxn_idx_tmp) < 0;
metlist = model.mets(met_idx_tmp);
slist = full(model.S(met_idx_tmp,rxn_idx_tmp));
flux_tmp = flux(ismember(model.rxns,'unmodeled_cofactor_formation'));

for i = 1:length(metlist)
    type_tmp = metlist(i);
    conc_tmp = -1 * slist(i) * flux_tmp / mu;
    if ~ismember(type_tmp,list)
        list = [list;type_tmp];
        concentration = [concentration;conc_tmp];
    else
        idx_tmptmp = ismember(list,type_tmp);
        concentration(idx_tmptmp) = concentration(idx_tmptmp) + conc_tmp;
    end
end


% 3) dummy cofactor

rxn_idx_tmp = ismember(model.rxns,'dilute_dummy_cofactor');
met_idx_tmp = model.S(:,rxn_idx_tmp) < 0;
metlist = model.mets(met_idx_tmp);
slist = full(model.S(met_idx_tmp,rxn_idx_tmp));
flux_tmp = flux(ismember(model.rxns,'dilute_dummy_cofactor'));

for i = 1:length(metlist)
    type_tmp = metlist(i);
    conc_tmp = -1 * slist(i) * flux_tmp / mu;
    if ~ismember(type_tmp,list)
        list = [list;type_tmp];
        concentration = [concentration;conc_tmp];
    else
        idx_tmptmp = ismember(list,type_tmp);
        concentration(idx_tmptmp) = concentration(idx_tmptmp) + conc_tmp;
    end
end









