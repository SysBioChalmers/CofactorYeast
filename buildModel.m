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
load('Protein_stoichiometry.mat');
model = addComplexFormationRxns(model,Protein_stoichiometry);

% manually update complex formation reactions for some complexes.
model = updateComplexformation(model);

%% Add complex dilution reactions
model = addComplexDilutionRxns(model);

%% Change the original biomass equation
% Estimate modeled proteome
f_modeled_protein = estimateModeledprotein(model);
f_modeled_protein = round(f_modeled_protein,2); %g/gProtein
% Estimate modeled cofactors
cofactor_info = estimateModeledcofactor(model);
save('cofactor_info.mat','cofactor_info');

% Add unmodeled cofactor reaction
model = addUnmodeledcofactor(model,cofactor_info,'unmodeled_cofactor[c]');

% Change the biomass equation
model = changeBiomass(model,f_modeled_protein,'r_4041','s_3717[c]','unmodeled_cofactor[c]','s_4206[c]','s_4205[c]');
% r_4041 is pseudo_biomass_rxn_id in the original GEM
% s_3717[c] is protein id in the original GEM
% unmodeled_cofactor[c] is newly added metabolite for unmodeled cofactor
% s_4206[c] is ion id in the original GEM
% s_4205[c] is cofactor id in the original GEM

%% Add dummy complex reactions
% Dummy complex is assumed to be a part of metabolic protein pool.
% Note that the dummy complex is synthesized or diluted in the unit of
% mmol/gCDW/h.

model = addDummy(model,'r_4047','s_3717[c]'); 
% r_4047 is pseudo_protein_rxn_id in the GEM
% s_3717[c] is protein id

model = addDummycofactor(model,cofactor_info);

%% Save model

save('CofactorYeast.mat','model');
writeCbModel(model,'xls','CofactorYeast.xls');

%% Collect kcats for enzymes

enzymedata = collectkcats(model);

% manually update kcats for some reactions.
enzymedata = updatekcats(enzymedata);
% CONFIDENCE SCORE 5 means maually assigned kcats, not completely correct.

%% Calculate molecular weight for each enzyme
enzymedata = calculateMW(enzymedata,ProteinSequence);

%% Collect cofactors for each enzyme
enzymedata = calculateCofactor(enzymedata,model);

save('enzymedata.mat','enzymedata');
toc;