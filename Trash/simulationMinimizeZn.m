%% simulationMinimizeZn
% Timing: ~  s
load('CofactorYeast.mat');
load('enzymedata.mat');

tic;

%% Set model
model = changeRxnBounds(model,'r_1714',-1000,'l'); %

% block some reactions
model = changeRxnBounds(model,'r_0886_1',0,'b'); % iso-reaction of PFK
model = changeRxnBounds(model,'r_4262_fwd',0,'b'); % citrate hydroxymutase
model = changeRxnBounds(model,'r_4262_rvs',0,'b'); % citrate hydroxymutase

% block some reactions that done in PMID: 28779005.
model = changeRxnBounds(model,'r_2045_rvs',0,'b'); % serine transport from [m] to [c]
model = changeRxnBounds(model,'r_0659_fwd',0,'b'); % isocitrate dehydrogenase (NADP)
model = changeRxnBounds(model,'r_0659_rvs',0,'b'); % isocitrate dehydrogenase (NADP)

% model = changeRxnBounds(model,'r_0725_fwd',0,'b'); % methenyltetrahydrofolate cyclohydrolase
% model = changeRxnBounds(model,'r_0918',0,'b'); % phosphoserine transaminase

model = changeRxnBounds(model,'r_4235',0,'b'); % weird reaction from glc to g6p

model = changeRxnBounds(model,'r_4216_rvs',0,'b'); % block the reaction to produce FMN without ATP

%% Add transport and exchange reactions
model = addReaction(model,'exchange_cd','reactionFormula','s_3783[e] -> ','reversible',true);
model = changeRxnBounds(model,'exchange_cd',-1000,'l');
model = addReaction(model,'transport_cd','reactionFormula','s_3783[e] -> s_3782[c]','reversible',true);



%% Set optimization
mu = 0.1;
model = changeRxnBounds(model,'r_2111',mu,'b');
tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
clear tot_protein f_modeled_protein;

%% Solve LPs
% Minimize glucose uptake as reference
% rxnID = 'r_1654'; % N source
rxnID = 'r_1714'; % C source
% rxnID = 'r_4596'; %minimize Zn uptake rate
osenseStr = 'Maximize';
model_tmp = model;
fileName = writeLP(model_tmp,mu,f,osenseStr,rxnID,enzymedata,1);
command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
system(command,'-echo');
[~,~,sol_full] = readSoplexResult('Simulation.lp.out',model_tmp);
glc = sol_full(strcmp(model_tmp.rxns,'r_1714'));
eth = sol_full(strcmp(model_tmp.rxns,'r_1761'));
co2 = sol_full(strcmp(model_tmp.rxns,'r_1672'));
o2 = sol_full(strcmp(model_tmp.rxns,'r_1992'));
[list, concentration] = simulatedCofactor(model_tmp,enzymedata,sol_full);

newzn = sol_full(strcmp(model_tmp.rxns,'r_1714'));
rxnID = 'r_4596'; 
osenseStr = 'Maximize';
model_tmp = model;
model_tmp = changeRxnBounds(model_tmp,'r_1714',newzn*1.000001,'b');
fileName = writeLP(model_tmp,mu,f,osenseStr,rxnID,enzymedata,1);
command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
system(command,'-echo');
[~,~,sol_full] = readSoplexResult('Simulation.lp.out',model_tmp);
glc = sol_full(strcmp(model_tmp.rxns,'r_1714'));
eth = sol_full(strcmp(model_tmp.rxns,'r_1761'));
co2 = sol_full(strcmp(model_tmp.rxns,'r_1672'));
o2 = sol_full(strcmp(model_tmp.rxns,'r_1992'));
[list, concentration] = simulatedCofactor(model_tmp,enzymedata,sol_full);

% Fix glucose uptake of reference and maximize dummy protein
% rxnID = 'dilute_dummy';
% osenseStr = 'Maximize';
% model_tmp = changeRxnBounds(model,'r_1714',-1.1-0.11,'l');
% model_tmp = changeRxnBounds(model_tmp,'r_1714',-1.1+0.11,'u');
% % model_tmp = changeRxnBounds(model_tmp,'r_1992',-3.1,'l');
% % model_tmp = changeRxnBounds(model_tmp,'r_1992',-2.5,'u');
% % model_tmp = changeRxnBounds(model_tmp,'r_1672',2.5,'l');
% % model_tmp = changeRxnBounds(model_tmp,'r_1672',3.1,'u');
% disp('Maximize dummy protein');
% fileName = writeLP(model_tmp,mu,f,osenseStr,rxnID,enzymedata,1);
% command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
% system(command,'-echo');
% [~,~,sol_full] = readSoplexResult('Simulation.lp.out',model_tmp);
% dummy_ref = sol_full(strcmp(model_tmp.rxns,'dilute_dummy'),:);

% % Minimize Zn uptake
% % rxnID = 'r_4596'; %minimize Zn uptake rate
% rxnID = 'r_1714'; %minimize glucose uptake rate
% osenseStr = 'Maximize';
% % model_tmp = changeRxnBounds(model,'dilute_dummy',0.0328,'b');
% model_tmp = changeRxnBounds(model,'dilute_dummy',0.0347,'b');
% disp('Minimize Zn uptake');
% fileName = writeLP(model_tmp,mu,f,osenseStr,rxnID,enzymedata,1);
% command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
% system(command,'-echo');
% [~,sol_status,sol_full] = readSoplexResult('Simulation.lp.out',model_tmp);
% 
% mu = sol_full(strcmp(model_tmp.rxns,'r_2111'),:);
% glc = -1*sol_full(strcmp(model_tmp.rxns,'r_1714'),:);
% etoh = sol_full(strcmp(model_tmp.rxns,'r_1761'),:);
% o2 = -1*sol_full(strcmp(model_tmp.rxns,'r_1992'),:);
% co2 = sol_full(strcmp(model_tmp.rxns,'r_1672'),:);
% % ac = fluxes(strcmp(model_tmp.rxns,'r_1634'),:);
% % aldh = fluxes(strcmp(model_tmp.rxns,'r_1631'),:);


toc;