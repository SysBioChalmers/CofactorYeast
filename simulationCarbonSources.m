%% simulationCS
% Timing: ~  s
load('CofactorYeast.mat');
load('enzymedata.mat');

tic;

%% Carbon sources
sCS_res = struct();
sCS_res.cslist = {'Glucose' 'Acetate' 'Ethanol' 'Fructose' 'Galactose' 'Glycerol' 'Maltose' 'Sucrose'};
exch_rxn_list = {'r_1714'  'r_1634'  'r_1761'  'r_1709'   'r_1710'    'r_1808'   'r_1931'  'r_2058'};

%% Set model

% block some reactions
model = changeRxnBounds(model,'r_0886_1',0,'b'); % iso-reaction of PFK
model = changeRxnBounds(model,'r_4262_fwd',0,'b'); % citrate hydroxymutase
model = changeRxnBounds(model,'r_4262_rvs',0,'b'); % citrate hydroxymutase

% block some reactions that done in PMID: 28779005.
model = changeRxnBounds(model,'r_2045_rvs',0,'b'); % serine transport from [m] to [c]
model = changeRxnBounds(model,'r_0659_fwd',0,'b'); % isocitrate dehydrogenase (NADP)
model = changeRxnBounds(model,'r_0659_rvs',0,'b'); % isocitrate dehydrogenase (NADP)

% model = changeRxnBounds(model,'r_0725_fwd',0,'b'); % methenyltetrahydrofolate cyclohydrolase
% model = changeRxnBounds(model,'r_0918',0,'b'); % phosphoserine transaminase

model = changeRxnBounds(model,'r_4235',0,'b'); % weird reaction from glc to g6p

model = changeRxnBounds(model,'r_4216_rvs',0,'b'); % block the reaction to produce FMN without ATP

%% Add transport and exchange reactions
model = addReaction(model,'exchange_cd','reactionFormula','s_3783[e] -> ','reversible',true);
model = changeRxnBounds(model,'exchange_cd',-1000,'l');
model = addReaction(model,'transport_cd','reactionFormula','s_3783[e] -> s_3782[c]','reversible',true);



%% Set optimization
rxnID = 'dilute_dummy';
osenseStr = 'Maximize';

tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
clear tot_protein f_modeled_protein;

%% Solve LPs
model = changeRxnBounds(model,'r_1714',0,'b'); % block default glucose uptake

sCS_res.mulist = zeros(1,length(sCS_res.cslist));
sCS_res.fluxes = zeros(length(model.rxns),length(sCS_res.cslist));

for i = 1:length(sCS_res.cslist)
    exrxn = exch_rxn_list{i};
    disp(['carbon source: ' sCS_res.cslist{i}]);
    model_tmp = changeRxnBounds(model,exrxn,-1000,'l');
    [mu_tmp,sol_full_tmp] = searchMaxgrowth(model_tmp,f,osenseStr,rxnID,enzymedata,1e-3);
    sCS_res.mulist(1,i) = mu_tmp;
    sCS_res.fluxes(:,i) = sol_full_tmp;
end

cd Results/;
save('sCS_res.mat','sCS_res');
cd ../;
clear;


%% Plot

toc;

