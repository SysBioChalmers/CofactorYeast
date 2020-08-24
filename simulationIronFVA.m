%% simulationIronFVA
% 
% Timing: ~ 1800 s
tic;
load('CofactorYeast.mat');
load('enzymedata.mat');

%% Set model
% set medium
model = setMedia(model,1);% minimal media (Delft media)
% set carbon source
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

factor_k_withoutcofator = 0;

%% Solve LPs
% reference
[mu_ref,flux_ref] = searchMaxgrowth(model,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator,1e-6);
q_fe_ref = flux_ref(strcmp(model.rxns,'r_1861'),1);

%% FVA
rxnID = 'r_1861';
mu = mu_ref;
model = changeRxnBounds(model,'r_2111',mu,'b');

osenseStr = 'Maximize';
disp(osenseStr);
fileName = writeLP(model,mu,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator);
command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
system(command,'-echo');
[q_fe_max,~,flux_max] = readSoplexResult('Simulation.lp.out',model);

osenseStr = 'Minimize';
disp(osenseStr);
fileName = writeLP(model,mu,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator);
command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
system(command,'-echo');
[q_fe_min,~,flux_min] = readSoplexResult('Simulation.lp.out',model);

sIFVA_res.flux_ref = flux_ref;
sIFVA_res.flux_max = flux_max;
sIFVA_res.flux_min = flux_min;

cd Results/;
save('sIFVA_res.mat','sIFVA_res');
cd ../;

toc;


