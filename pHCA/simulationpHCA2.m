%% simulationpHCA2

% Timing: ~ 50000 s
tic;

% load model
load('modelpHCA.mat');
load('enzymedatapHCA.mat');

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
rxnID = 'new_r_HCAex'; %maximize noscapine production rate
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
load('spHCA1_fluxes.mat');
q_fe_tmp = fluxes(strcmp(model.rxns,'r_1861'),:);
q_fe_tmp = round(min(q_fe_tmp),6);
clear fluxes;

q_fe_list = q_fe_tmp*(0.02:0.02:1);
clear q_fe_tmp;

init_mu = 0;
step_mu = 0.01;
fluxes = zeros(length(model.rxns),0);
record = zeros(1,0);
for i = 1:length(q_fe_list)
    q_fe = q_fe_list(i);
    model_tmp = model;
    model_tmp.lb(strcmp(model.rxns,'r_1861')) = q_fe;
    
    sol_status = 'optimal';
    mu = init_mu;
    while strcmp(sol_status,'optimal')
        mu = mu + step_mu;
        
        if mu == 0.38
            mu = 0.379;
        end
        
        model_tmp_tmp = changeRxnBounds(model_tmp,'r_2111',mu,'b');
        disp(['q_Fe = ' num2str(q_fe) '; mu = ' num2str(mu)]);
        fileName = writeLP(model_tmp_tmp,mu,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator);
        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        [~,sol_status,sol_full] = readSoplexResult('Simulation.lp.out',model_tmp_tmp);
        if strcmp(sol_status,'optimal')
            fluxes = [fluxes sol_full];
            record = [record q_fe];
        end
    end
end
spHCA2_res.fluxes = fluxes;
spHCA2_res.record = record;

save('spHCA2_res.mat','spHCA2_res');

clear;

toc;

