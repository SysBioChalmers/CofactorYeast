%% Simulation

load('CofactorYeast.mat');
load('enzymedata.mat');

%% Set model


sol = optimizeCbModel(model);


%% Set optimization
rxnID = 'r_1714'; %minimize glucose uptake rate
osenseStr = 'Maximize';

tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id
f = tot_protein * f_modeled_protein;
clear tot_protein f_modeled_protein;

mu = 0.01;
fileName = writeLP(model,mu,f,osenseStr,rxnID,enzymedata,1);

command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
system(command,'-echo');
