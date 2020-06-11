load('sB_res.mat');

load('CofactorYeast.mat');

ion_id_list = {'CA';'CU';'FE';'K';'MG';'MN';'NA';'ZN'};
ion_ex_list = {'r_4600'; ... % Ca(2+) exchange
               'r_4594'; ... % Cu(2+) exchange
               'r_1861'; ... % iron(2+) exchange
               'r_2020'; ... % potassium exchange
               'r_4597'; ... % Mg(2+) exchange
               'r_4595'; ... % Mn(2+) exchange
               'r_2049'; ... % sodium exchange
               'r_4596'};... % Zn(2+) exchange

idx_gpr = ~ismember(model.grRules,'');
idx_zero_flux = any(sB_res.fluxes,2);
idx_combined = idx_gpr & idx_zero_flux;

r_values = zeros(678,1);
p_values = zeros(678,1);

for i = 1:length(ion_id_list)
    idx_ion_ex_tmp = ismember(model.rxns,ion_ex_list(i));
    ion_uptake_tmp = -1*sB_res.fluxes(idx_ion_ex_tmp,:);
    non_zero_fluxes = sB_res.fluxes(idx_combined,:);
    [a, ~] = size(non_zero_fluxes);
    for j = 1:a
        [r_tmp,p_tmp] = corrcoef(ion_uptake_tmp,non_zero_fluxes(j,:));
        r_values(j) = r_tmp(1,2);
        p_values(j) = p_tmp(1,2);
        scatter(ion_uptake_tmp,non_zero_fluxes(j,:));
    end
    
end

