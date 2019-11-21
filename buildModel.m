%% buildModel

%% Import yeast 7.6

% load('Yeast7.6.mat');

%% Import yeast8
% cd GEMs/yeast_8.3.4/;
% org_model = readCbModel('yeastGEM.xml');
% cd ../../;
% save('Yeast8.mat','org_model');
load('Yeast8.mat');

%% Modify the original model
model_updated = org_model;
% after updated the model, protein list in the file
% "ProteinStoichiometry.xlsx" should also be updated, and
% "process_protein_stoichiometry.m" should be re-run.
clear org_model;

%% Reformulate the orginal model
%   1. Reactions with isozymes (i.e., 'OR' case in GPR) will be copied,
%   each isozyme will be assigned to each copy. 
%   2. Reversible reactions will be split into forward and reverse ones.

model_split = splitModel(model_updated);
save('model_split.mat','model_split');
clear model_updated;

%% Add translation reactions
load('ProteinSequence.mat');
model = addTranslationRxns(model_split,ProteinSequence);

%% Add cofactor binding reactions
model = addCofactorRxns(model);

%% Add complex formation reactions based on protein stoichiometry
% promiscuous = findPromiscuous(model_split);
load('Determined_stoichiometry.mat');
model = addComplexFormationRxns(model,Determined_stoichiometry);

% manually update complex formation reactions for some complexes.
% model = updateComplexformation(model);


%% Collect kcats for complexes

enzymedata = collectkcats(model);

% manually update kcats for some reactions.
% enzymedata = updatekcats(enzymedata);


