%% simulationIron
% Timing: ~ 9500 s
tic;
load('CofactorYeastExpand.mat');
load('enzymedataExpand.mat');

%% Set model
% set medium
model = setMedia(model,1);% minimal media (Delft media)
% set carbon source
model = changeRxnBounds(model,'r_1714',-1000,'l');% glucose
% set oxygen
model = changeRxnBounds(model,'r_1992',-1000,'l');
% block reactions
model = blockRxns(model);
% model = changeRxnBounds(model,'r_1631',0,'b');% acetaldehyde production

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
% reference
[~,flux_ref] = searchMaxgrowth4Expand(model,f,osenseStr,rxnID,enzymedata,1e-6,1,0);
q_fe_ref = flux_ref(strcmp(model.rxns,'r_1861'),1);

%% Solve LPs
k_cf = [0 0.1 0.2 0.5 0.9];
lower_fe = 0.8;

fluxes = zeros(length(model.rxns),length(k_cf));
for i = 1:length(k_cf)
    disp([num2str(i) '/' num2str(length(k_cf))]);
    model_tmp = changeRxnBounds(model,'r_1861',q_fe_ref*lower_fe,'b');
    [~,flux_tmp] = searchMaxgrowth4Expand(model_tmp,f,osenseStr,rxnID,enzymedata,1e-6,1,k_cf(i));
    fluxes(:,i) = flux_tmp;
end

sI_res.k_cf = k_cf;
sI_res.fluxes = fluxes;
sI_res.flux_ref = flux_ref;

cd Results/;
save('sI_res.mat','sI_res');
cd ../;
clear;

toc;

