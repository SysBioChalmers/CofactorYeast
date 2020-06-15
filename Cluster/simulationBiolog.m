%% simulationBiolog
function simulationBiolog(a,b)
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
% Upper bound of source uptake, estimated by reference condition.
N_ub = 2.6*0.6/0.4;
P_ub = 0.11*0.6/0.4;
S_ub = 0.042*0.6/0.4;

[~,raw,~] = xlsread('Exp_biolog.xlsx');
metNameID = raw(2:end,2);
sourceType = raw(2:end,3);
expRes = raw(2:end,4);
simRes = raw(2:end,5);
idxG = ismember(expRes,'G') & ismember(simRes,'G');

metNameID = metNameID(idxG);
sourceType = sourceType(idxG);

fluxes = zeros(length(model.rxns),b-a+1);
labels = cell(1,b-a+1);

for i = a:b
    metNameID_tmp = metNameID{i};
    metE = strcat(metNameID_tmp,' [extracellular]');
    idx_E = model.S(ismember(model.metNames,metE),:)' ~= 0 & findExcRxns(model);
    exRxnID = model.rxns(idx_E);
    metID_i = model.mets(ismember(model.metNames,metE));
    metID_tmp = strrep(metID_i,'[','_');
    metID_tmp = strrep(metID_tmp,']','');
    source_tmp = sourceType{i};
    label_tmp = strcat(source_tmp,'_',metID_tmp);
    
    [Ematrix, elements] = computeElementalMatrix(model,metID_i,true,true);
    if strcmp(source_tmp,'C')
        model_tmp = changeRxnBounds(model,'r_1714',0,'l');% glucose
        model_tmp = changeRxnBounds(model_tmp,cell2mat(exRxnID),-1000,'l');
    elseif strcmp(source_tmp,'N')
        model_tmp = changeRxnBounds(model,'r_1654',0,'l');% ammonium
        bound_tmp = Ematrix(ismember(elements,{'N'})) * N_ub;
        model_tmp = changeRxnBounds(model_tmp,cell2mat(exRxnID),-1*bound_tmp,'l');
    elseif strcmp(source_tmp,'P')
        model_tmp = changeRxnBounds(model,'r_2005',0,'l');% phosphate
        bound_tmp = Ematrix(ismember(elements,{'P'})) * P_ub;
        model_tmp = changeRxnBounds(model_tmp,cell2mat(exRxnID),-1*bound_tmp,'l');
    elseif strcmp(source_tmp,'S')
        model_tmp = changeRxnBounds(model,'r_2060',0,'l');% sulphate
        bound_tmp = Ematrix(ismember(elements,{'S'})) * S_ub;
        model_tmp = changeRxnBounds(model_tmp,cell2mat(exRxnID),-1*bound_tmp,'l');
    end
    
    [~,flux_tmp] = searchMaxgrowthSpecial(model_tmp,cell2mat(label_tmp),f,osenseStr,rxnID,enzymedata,1e-6);
    fluxes(:,i-a+1) = flux_tmp;
    labels(1,i-a+1) = label_tmp;
end

filename_1 = strcat('sB_fluxes_',num2str(a),'.mat');
filename_2 = strcat('sB_labels_',num2str(a),'.mat');
cd tmp_results/;
save(filename_1,'fluxes');
save(filename_2,'labels');
cd ../;
end
