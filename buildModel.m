%% buildModel
% Timing: ~ 500 s

tic;
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

%% Add complex dilution reactions
model = addComplexDilutionRxns(model);

%% Add dummy complex reactions
% Dummy complex is assumed to be a part of metabolic protein pool.
% Note that the dummy complex is synthesized or diluted in the unit of
% g/gCDW/h.

model = addDummy(model,'r_4047','s_3717[c]'); 
% r_4047 is pseudo_protein_rxn_id in the GEM
% s_3717[c] is protein id


%% Change protein in the biomass equation
% Estimate modeled proteome
f_modeled_protein = estimateModeledprotein(model);
f_modeled_protein = round(f_modeled_protein,2); %g/gProtein

% Change the coefficient of the protein in the biomass equation to
% unmodeled fraction
model = changeBiomass(model,f_modeled_protein,'r_4041','s_3717[c]');
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

save('CofactorYeast.mat','model');
writeCbModel(model,'xls','CofactorYeast.xls');

%% Collect kcats for enzymes

enzymedata = collectkcats(model);

% manually update kcats for some reactions.
% enzymedata = updatekcats(enzymedata);
% CONFIDENCE SCORE 4 means maually collected kcats.

%% Calculate molecular weight for each enzyme
enzymedata = calculateMW(enzymedata,ProteinSequence);

save('enzymedata.mat','enzymedata');
toc;