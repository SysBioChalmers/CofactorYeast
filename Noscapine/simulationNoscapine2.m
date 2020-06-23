%% simulationNoscapine2

% Timing: ~ 16000 s
tic;

load('modelNoscapine.mat');
load('enzymedataNoscapine.mat');

%% Set model
% set medium
model = setMedia(model,2);% yeast nitrogen base without amino acids
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
q_fe_list = -3.5e-4*(0.05:0.05:1);
factor_k_withoutcofator = 0.2;
init_mu = 0;
step_mu = 0.02;
fluxes = zeros(length(model.rxns),0);
record = zeros(1,0);
for i = 1:length(q_fe_list)
    q_fe = q_fe_list(i);
    model_tmp = model;
    model_tmp.lb(ismember(model_tmp.rxnNames,'iron(3+) exchange')) = q_fe;
    
    sol_status = 'optimal';
    mu = init_mu;
    while strcmp(sol_status,'optimal')
        mu = mu + step_mu;
        model_tmp_tmp = changeRxnBounds(model_tmp,'r_2111',mu,'b');
        disp(['q_Fe = ' num2str(q_fe) '; mu = ' num2str(mu)]);
        fileName = writeLP4Expand(model_tmp_tmp,mu,f,osenseStr,rxnID,enzymedata,1,factor_k_withoutcofator);
        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        [~,sol_status,sol_full] = readSoplexResult('Simulation.lp.out',model_tmp_tmp);
        if strcmp(sol_status,'optimal')
            fluxes = [fluxes sol_full];
            record = [record q_fe];
        end
    end
end
sN2_res.fluxes = fluxes;
sN2_res.record = record;

cd Results/;
save('sN2_res.mat','sN2_res');
cd ../;
clear;

toc;

