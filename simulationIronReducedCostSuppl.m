%% simulationIronReducedCostSuppl
% Add palmitoleate, oleate and ergosterol to see how growth changes

% Timing: ~ 15000 s
tic;
load('CofactorYeast.mat');
load('enzymedata.mat');

%% Set model
% set medium
model = setMedia(model,1);% minimal media (Delft media)
% set carbon source
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

%% Solve LPs
% reference
[~,flux_ref100] = searchMaxgrowth(model,f,f_mito,osenseStr,rxnID,enzymedata,0,1e-6);
q_fe_ref100 = flux_ref100(strcmp(model.rxns,'r_1861'),1);

%% Solve LPs

factor_k_withoutcofator = 0.5;
lower_fe = 0.5;
q_fe_new = q_fe_ref100*lower_fe;

model_fe50 = changeRxnBounds(model,'r_1861',q_fe_new,'b');
[~,flux_ref50] = searchMaxgrowth(model_fe50,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator,1e-6);


met_Exchanges = {'r_1994' ... % palmitoleate exchange
                 'r_2189' ... % oleate exchange
                 'r_1757'};   %	ergosterol exchange
met_IDs =   {'Palmitoleate' 'Oleate' 'Ergosterol'};
q_delta = -0.001;

fluxes = zeros(length(model.rxns),length(met_Exchanges));
for i = 1:length(met_Exchanges)
    model_tmp = model_fe50;
    model_tmp = changeRxnBounds(model_tmp,met_Exchanges{i},q_delta,'l');
    disp(['Adding ' met_IDs{i}]);
    [~,flux_tmp] = searchMaxgrowth(model_tmp,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator,1e-6);
    fluxes(:,i) = flux_tmp;
end

sIRCS_res.lables = ['Ref' met_IDs];
sIRCS_res.fluxes = [flux_ref50 fluxes];

cd Results/;
save('sIRCS_res.mat','sIRCS_res');
cd ../;
clear;


%% simulationRefReducedCostAA
% Add each amino acid to see how growth changes

load('CofactorYeast.mat');
load('enzymedata.mat');

%% Set model
% set medium
model = setMedia(model,1);% minimal media (Delft media)
% set carbon source
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

%% Solve LPs
% reference
[~,flux_ref] = searchMaxgrowth(model,f,f_mito,osenseStr,rxnID,enzymedata,0,1e-6);

%% Solve LPs

factor_k_withoutcofator = 0;

met_Exchanges = {'r_1994' ... % palmitoleate exchange
                 'r_2189' ... % oleate exchange
                 'r_1757'};   %	ergosterol exchange
met_IDs =   {'Palmitoleate' 'Oleate' 'Ergosterol'};
q_delta = -0.001;

fluxes = zeros(length(model.rxns),length(met_Exchanges));
for i = 1:length(met_Exchanges)
    model_tmp = model;
    model_tmp = changeRxnBounds(model_tmp,met_Exchanges{i},q_delta,'l');
    disp(['Adding ' met_IDs{i}]);
    [~,flux_tmp] = searchMaxgrowth(model_tmp,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator,1e-6);
    fluxes(:,i) = flux_tmp;
end

sRRCS_res.lables = ['Ref' met_IDs];
sRRCS_res.fluxes = [flux_ref fluxes];

cd Results/;
save('sRRCS_res.mat','sRRCS_res');
cd ../;
clear;

toc;