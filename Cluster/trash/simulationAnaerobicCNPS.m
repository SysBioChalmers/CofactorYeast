%% simulationAnaerobicCNPS
function simulationAnaerobicCNPS(a,b)
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

%% Set anaerobic condition (according to yeast8 repo)
%1st change: Refit GAM and NGAM to exp. data
GAM   = 30.49;  %Data from Nissen et al. 1997
NGAM  = 0;      %Refit done in Jouthen et al. 2012
bioPos = strcmp(model.rxnNames,'biomass pseudoreaction');
for i = 1:length(model.mets)
    S_ix  = model.S(i,bioPos);
    isGAM = sum(strcmp({'ATP [cytoplasm]','ADP [cytoplasm]','H2O [cytoplasm]', ...
        'H+ [cytoplasm]','phosphate [cytoplasm]'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        model.S(i,bioPos) = sign(S_ix)*GAM;
    end
end
model = changeRxnBounds(model,'r_4046',NGAM,'l');

%2nd change: Removes the requirement of heme a, NAD(PH), coenzyme A in the biomass equation
%            (not used under anaerobic conditions)
mets = {'s_3714[c]','s_1198[c]','s_1203[c]','s_1207[c]','s_1212[c]','s_0529[c]'};
[~,met_index] = ismember(mets,model.mets);
model.S(met_index,strcmp(model.rxns,'r_4598')) = 0;

%3rd change: Changes media to anaerobic (no O2 uptake and allows sterol
%            and fatty acid exchanges)
model.lb(strcmp(model.rxns,'r_1992')) = 0;        %O2
model.lb(strcmp(model.rxns,'r_1757')) = -1000;    %ergosterol
model.lb(strcmp(model.rxns,'r_1915')) = -1000;    %lanosterol
model.lb(strcmp(model.rxns,'r_1994')) = -1000;    %palmitoleate
model.lb(strcmp(model.rxns,'r_2106')) = -1000;    %zymosterol
model.lb(strcmp(model.rxns,'r_2134')) = -1000;    %14-demethyllanosterol
model.lb(strcmp(model.rxns,'r_2137')) = -1000;    %ergosta-5,7,22,24(28)-tetraen-3beta-ol
model.lb(strcmp(model.rxns,'r_2189')) = -1000;    %oleate

%4th change: Blocked pathways for proper glycerol production
%Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
model.ub(strcmp(model.rxns,'r_0713_rvs')) = 0; %Mithocondria
model.ub(strcmp(model.rxns,'r_0714_rvs')) = 0; %Cytoplasm
%Block glycerol dehydroginase (only acts in microaerobic conditions)
model.ub(strcmp(model.rxns,'r_0487')) = 0;
%Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
model.ub(strcmp(model.rxns,'r_0472')) = 0;

%5th change: Set free some reactions that are blocked in aerobic conditions
model = changeRxnBounds(model,'r_0487',1000,'u');
model.ub(ismember(model.rxns,'r_0487_withoutcofactor')) = 1000;

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
idxG = ~ismember(metNameID,'');
metNameID = metNameID(idxG);
sourceType = sourceType(idxG);

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
    label_tmp = strcat(source_tmp,'_',metID_tmp,'_Anaerobic');
    
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

filename_1 = strcat('sAnaerobicCNPS_fluxes_',num2str(a),'.mat');
filename_2 = strcat('sAnaerobicCNPS_labels_',num2str(a),'.mat');
cd tmp_results/;
save(filename_1,'fluxes');
save(filename_2,'labels');
cd ../;
end
