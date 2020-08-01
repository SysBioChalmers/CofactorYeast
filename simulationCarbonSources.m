%% simulationCS
% Timing: ~ 9000 s
load('CofactorYeast.mat');
load('enzymedata.mat');
tic;

%% Carbon sources
sCS_res = struct();
sCS_res.cslist = {'Glucose' 'Acetate' 'Ethanol' 'Fructose' 'Galactose' 'Glycerol' 'Maltose' 'Sucrose'};
exch_rxn_list = {'r_1714'  'r_1634'  'r_1761'  'r_1709'   'r_1710'    'r_1808'   'r_1931'  'r_2058'};

%% Set model
% set medium
model = setMedia(model,2);% yeast nitrogen base without amino acids
% set oxygen
model = changeRxnBounds(model,'r_1992',-1000,'l');
% block reactions
model = blockRxns(model);
model = changeRxnBounds(model,'r_1631',0,'b');% acetaldehyde production

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

sCS_res.mulist = zeros(1,length(sCS_res.cslist));
sCS_res.fluxes = zeros(length(model.rxns),length(sCS_res.cslist));

for i = 1:length(sCS_res.cslist)
    exrxn = exch_rxn_list{i};
    disp(['carbon source: ' sCS_res.cslist{i}]);
    model_tmp = changeRxnBounds(model,exrxn,-1000,'l');
    [mu_tmp,sol_full_tmp] = searchMaxgrowth(model_tmp,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator,1e-6);
    sCS_res.mulist(1,i) = mu_tmp;
    sCS_res.fluxes(:,i) = sol_full_tmp;
end

cd Results/;
save('sCS_res.mat','sCS_res');
cd ../;
clear;

toc;

