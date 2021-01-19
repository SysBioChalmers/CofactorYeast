%% construct pHCA pathway
load('CofactorYeast.mat');
load('enzymedata.mat');
[model,enzymedata] = addNewPathway(model,enzymedata,'New_pathway_pHCA.xlsx');
[model,enzymedata] = expandModel(model,enzymedata);
save('modelpHCA.mat','model');
save('enzymedatapHCA.mat','enzymedata');