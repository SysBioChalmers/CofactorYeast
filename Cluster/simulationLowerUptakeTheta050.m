%% simulationLowerUptakeTheta050
function simulationLowerUptakeTheta050(ion_mw)
addpath(genpath('../../CofactorYeast/')); %add to path
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/libSBML-5.15.0-matlab'))
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools'))
workdir = pwd;
cd '/cephyr/users/feiranl/Hebbe/tools/cobratoolbox'
initCobraToolbox
savepath '~/pathdef.m'
cd(workdir)

load('CofactorYeast.mat');
load('enzymedata.mat');

%% Set model
% set medium
model = setMedia(model,1);% minimal media (Delft media)
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

factor_k_withoutcofator = 0.5;

%% Solve LPs
ion_mw_list = [39 24 56 65 40 55 64 23];
ion_list = {'K' 'MG' 'FE' 'ZN' 'CA' 'MN' 'CU' 'NA'};
rxnName_list = {'potassium exchange'; ...
                'Mg(2+) exchange'; ...
                'iron(2+) exchange'; ...
                'Zn(2+) exchange'; ...
                'Ca(2+) exchange'; ...
                'Mn(2+) exchange'; ...
                'Cu2(+) exchange'; ...
                'sodium exchange'};
decrease_value = 0.95:-0.05:0.05;

i = ion_mw_list == ion_mw;
ion = ion_list{i};

fluxes = zeros(length(model.rxns),length(decrease_value)+1);
labels = cell(1,length(decrease_value)+1);

% reference
[~,fluxes_ref] = searchMaxgrowthSpecial(model,strcat(ion,'_1'),f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator,1e-6);

fluxes(:,1) = fluxes_ref;
labels(1,1) = {strcat(ion,'_1')};

% decreased uptake

rxn_tmp = rxnName_list(i);
ref_uptake = fluxes_ref(ismember(model.rxnNames,rxn_tmp));
for j = 1:length(decrease_value)
    new_uptake = ref_uptake*decrease_value(j);
    model_tmp = model;
    model_tmp.lb(ismember(model.rxnNames,rxn_tmp)) = new_uptake;
    model_tmp.ub(ismember(model.rxnNames,rxn_tmp)) = new_uptake;
    str_tmp = num2str(decrease_value(j));
    str_tmp = strrep(str_tmp,'.','_');
    label_tmp = strcat(ion,'_',str_tmp);
    [~,fluxes_tmp] = searchMaxgrowthSpecial(model_tmp,label_tmp,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator,1e-6);
    fluxes(:,j+1) = fluxes_tmp;
    labels(1,j+1) = {label_tmp};
end

filename_1 = strcat('sLUT050_fluxes_',ion,'.mat');
filename_2 = strcat('sLUT050_labels_',ion,'.mat');
cd tmp_results/;
save(filename_1,'fluxes');
save(filename_2,'labels');
cd ../;
end