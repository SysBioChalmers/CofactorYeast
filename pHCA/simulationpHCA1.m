%% simulationpHCA1

% Timing: ~ 850 s
tic;

% load model
load('modelpHCA.mat');
load('enzymedatapHCA.mat');

soplexpath = '/Users/cheyu/build/bin/soplex'; % change this to the soplex path on your PC


%% Set model
% set medium
model = setMedia(model,1);% minimal medium
% set carbon source
model = changeRxnBounds(model,'r_1714',-1000,'l');% glucose
% set oxygen
model = changeRxnBounds(model,'r_1992',-1000,'l');
% block reactions
model = blockRxns(model);


%% Set optimization
rxnID = 'new_r_HCAex'; %maximize pHCA production rate
osenseStr = 'Maximize';

tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
f_mito = 0.1;
clear tot_protein f_modeled_protein;

factor_k_withoutcofator = 0.5;

%% Solve LPs
mu_list = [0.02:0.02:0.36 0.379];

fluxes = zeros(length(model.rxns),length(mu_list));

for i = 1:length(mu_list)
    mu = mu_list(i);
    model_tmp = changeRxnBounds(model,'r_2111',mu,'b');
    disp(['mu = ' num2str(mu)]);
    fileName = writeLP(model_tmp,mu,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator);
    command = sprintf([soplexpath,' -s0 -g5 -t300 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s'],fileName,fileName);
    system(command,'-echo');
    [sol_obj,sol_status,sol_full] = readSoplexResult('Simulation.lp.out',model_tmp);
    disp(['solution status: ' sol_status]);
    if strcmp(sol_status,'optimal')
        fluxes(:,i) = sol_full;
    end
end

save('spHCA1_fluxes.mat','fluxes');
clear;

toc;

