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

%% Formulate S matrix for complex formation based on protein stoichiometry
% promiscuous = findPromiscuous(model_split);




