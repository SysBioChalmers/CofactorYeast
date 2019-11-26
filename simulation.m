%% Simulation

load('CofactorYeast.mat');
load('enzymedata.mat');
% enzymedata = updatekcats(enzymedata);
tic;

%% Set model

model = changeRxnBounds(model,'r_1714',-1000,'l'); % 
model = changeRxnBounds(model,'r_0886_1',0,'b'); % block iso-reaction of PFK

%% Set optimization
rxnID = 'r_1714'; %minimize glucose uptake rate
% rxnID = 'dilute_dummy'; %minimize glucose uptake rate
osenseStr = 'Maximize';

tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
clear tot_protein f_modeled_protein;

mu = 0.2;
model = changeRxnBounds(model,'r_2111',mu,'b');
fileName = writeLP(model,mu,f,osenseStr,rxnID,enzymedata,1);

%% Solve LP
command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
system(command,'-echo');
[sol_obj,sol_status,sol_full] = readSoplexResult('Simulation.lp.out',model);

% dummy complex
if strcmp(sol_status,'optimal')
    dummy = sol_full(strcmp(model.rxns,'dilute_dummy'))/mu; %g/gCDW
    etoh = sol_full(strcmp(model.rxns,'r_1761'));
    glc = -1*sol_full(strcmp(model.rxns,'r_1714'));
end

toc;

