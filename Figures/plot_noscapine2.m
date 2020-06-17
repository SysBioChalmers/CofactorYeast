load('sN2_res.mat');
load('modelNoscapine.mat');

q_fe_list = unique(sN2_res.record);
tot_raw = sum(sN2_res.record == min(q_fe_list));
data_mu = zeros(tot_raw,length(q_fe_list));
data_Noscapine = zeros(tot_raw,length(q_fe_list));

for i = 1:length(q_fe_list)
    idx_tmp = sN2_res.record == q_fe_list(i);
    flux_tmp = sN2_res.fluxes(:,idx_tmp);
    mu = flux_tmp(strcmp(model.rxns,'r_2111'),:);
    glc = -1*flux_tmp(strcmp(model.rxns,'r_1714'),:);
    Noscapine = flux_tmp(strcmp(model.rxns,'new_r_eNoscapine'),:);
    Noscapine = Noscapine./glc;
    data_mu(:,i) = [mu nan(1,tot_raw-length(mu))];
    data_Noscapine(:,i) = [Noscapine nan(1,tot_raw-length(Noscapine))];
end

surf(-q_fe_list,data_mu,data_Noscapine)
