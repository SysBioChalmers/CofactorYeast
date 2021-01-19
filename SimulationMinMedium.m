%% simulationMinMedium
% Timing: ~ 1600 s
load('CofactorYeast.mat');
load('enzymedata.mat');
tic;

%% Set model
% set medium
model = setMedia(model,1);% minimal media (Delft media) (default)
model = changeRxnBounds(model,'r_1714',-1000,'l');% glucose
% set oxygen
model = changeRxnBounds(model,'r_1992',-1000,'l');
% block reactions
model = blockRxns(model);

%% Set optimization
rxnID = 'dilute_dummy';
osenseStr = 'Maximize';

tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
f_mito = 0.1;
clear tot_protein f_modeled_protein;

factor_k_withoutcofator = 0;

%% Solve LPs

sMM_res = struct();
sMM_res.mulist = 0;
sMM_res.fluxes = zeros(length(model.rxns),1);

[mu_tmp,sol_full_tmp] = searchMaxgrowth(model,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator,1e-6);
sMM_res.mulist(1,1) = mu_tmp;
sMM_res.fluxes(:,1) = sol_full_tmp;

cd Results/;
save('sMM_res.mat','sMM_res');
cd ../;
clear;

toc;

