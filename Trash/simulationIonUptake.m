%% simulationIonUptake
% Timing: ~ 65000 s
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
model = changeRxnBounds(model,'r_1631',0,'b');% acetaldehyde production ???????

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
ion_list = {'K' 'MG' 'FE' 'ZN' 'CA' 'MN' 'CU' 'NA'};
rxnName_list = {'potassium exchange'; ...
                'Mg(2+) exchange'; ...
                'iron(2+) exchange'; ...
                'Zn(2+) exchange'; ...
                'Ca(2+) exchange'; ...
                'Mn(2+) exchange'; ...
                'Cu2(+) exchange'; ...
                'sodium exchange'};

% reference
disp('Simulating reference state... ');
[~,fluxes_ref] = searchMaxgrowth(model,f,osenseStr,rxnID,enzymedata,1e-6);
sIU_res.fluxes_ref = fluxes_ref;

% decreased ion uptake
decrease_value = 0.8:-0.2:0;
sIU_res.fluxes = zeros(length(model.rxns),length(rxnName_list)*length(decrease_value));
sIU_res.ionlist = ion_list;

for i = 1:length(rxnName_list)
    rxn_tmp = rxnName_list(i);
    ref_uptake = fluxes_ref(ismember(model.rxnNames,rxn_tmp));
    for j = 1:length(decrease_value)
        disp(['Simulating ' ion_list{i} ': ' num2str(decrease_value(j))]);
        new_uptake = ref_uptake*decrease_value(j);
        model_tmp = model;
        model_tmp.lb(ismember(model.rxnNames,rxn_tmp)) = new_uptake;
        model_tmp.ub(ismember(model.rxnNames,rxn_tmp)) = new_uptake;
        [~,fluxes_tmp] = searchMaxgrowth(model_tmp,f,osenseStr,rxnID,enzymedata,1e-6);
        sIU_res.fluxes(:,(i-1)*length(decrease_value)+j) = fluxes_tmp;
    end
end

cd Results/;
save('sIU_res.mat','sIU_res');
cd ../;
clear;

toc;
