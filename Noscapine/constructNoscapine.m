%% construct noscapine pathway
load('CofactorYeast.mat');
load('enzymedata.mat');
[model,enzymedata] = addNewPathway(model,enzymedata,'New_pathway_noscapine.xlsx');
[model,enzymedata] = expandModel(model,enzymedata);%%%?????
save('modelNoscapine.mat','model');
save('enzymedataNoscapine.mat','enzymedata');