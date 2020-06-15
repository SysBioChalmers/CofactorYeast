%% simulationNoscapine

% Timing: ~ 700 s
tic;

% load original model
load('CofactorYeast.mat');
load('enzymedata.mat');

%% add new pathway
[model,enzymedata] = addNewPathway(model,enzymedata,'New_pathway_noscapine.xlsx');

%% Set model
% set medium
model = setMedia(model,1);% minimal media (Delft media)
% set carbon source
model = changeRxnBounds(model,'r_1714',-1000,'l');% glucose
% set oxygen
model = changeRxnBounds(model,'r_1992',-1000,'l');
% block reactions
model = blockRxns(model);
% model = changeRxnBounds(model,'r_1631',0,'b');% acetaldehyde production

%% Set optimization
rxnID = 'new_r_eNoscapine'; %maximize noscapine production rate
osenseStr = 'Maximize';

tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
clear tot_protein f_modeled_protein;

%% Solve LPs
% mu_list = 0.1:0.1:0.4;
mu_list = 0.01;
fluxes = zeros(length(model.rxns),length(mu_list));

for i = 1:length(mu_list)
    mu = mu_list(i);
    model_tmp = changeRxnBounds(model,'r_2111',mu,'b');
    disp(['mu = ' num2str(mu)]);
    fileName = writeLP(model_tmp,mu,f,osenseStr,rxnID,enzymedata,1);
    command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
    system(command,'-echo');
    [sol_obj,sol_status,sol_full] = readSoplexResult('Simulation.lp.out',model_tmp);
    disp(['solution status: ' sol_status]);
    if strcmp(sol_status,'optimal')
        fluxes(:,i) = sol_full;
    end
end


toc;


