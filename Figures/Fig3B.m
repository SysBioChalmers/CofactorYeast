%% Plot lower uptake of each metal ion
% PCA
load('sLU_res.mat');
load('sLUOCNPS_res.mat');
load('CofactorYeast.mat');

fluxes = [sLU_res.fluxes sLUOCNPS_res.fluxes];
labels = [sLU_res.labels sLUOCNPS_res.labels];
rxns = model.rxns;
clear sLU_res sLUOCNPS_res;
% remove non-metabolic fluxes
exclude = {'_translated','_cofactorbound','_enzyme_formation','_enzyme_dilution'};
metrxn = ~contains(model.rxns,exclude);
fluxes = fluxes(metrxn,:); 
rxns = rxns(metrxn,:);
clear exclude metrxn;
% remove the reactions with at least one flux > 100 or < -100
idx_over100 = any(abs(fluxes) > 100,2);
fluxes = fluxes(~idx_over100,:);
rxns = rxns(~idx_over100);
clear idx_over100;
% remove the reactions with zero flux at any simulation
idx_zero = all(fluxes == 0,2);
fluxes = fluxes(~idx_zero,:);
rxns = rxns(~idx_zero);
clear idx_zero;
% remove the reactions with small flux (<0.0001) at any simulation
idx_small = all(abs(fluxes) < 0.0001,2);
fluxes = fluxes(~idx_small,:);
rxns = rxns(~idx_small);
clear idx_small;

%% PCA plot
ion_id_list = {'CA';'CU';'FE';'K';'MG';'MN';'NA';'ZN'};
ion_clr_list = [228,26,28           % CA
                55,126,184          % CU
                255,127,0           % FE
                77,175,74           % K
                152,78,163          % MG
                166,86,40           % MN
                247,129,191         % NA
                153,153,153]/255;   % ZN
ref_clr = [0,0,0];

figure('Name','PCA');