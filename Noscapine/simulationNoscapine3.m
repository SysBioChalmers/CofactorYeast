%% simulationNoscapine3

% Set growth rate at 0.1/h, unlimited iron uptake and 50% ref iron uptake
% Set uptake rate of each amino acid at 0.01 then maximize noscapine.

% Timing: ~ 2800 s
tic;

% load model
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

%% Set optimization
rxnID = 'new_r_eNoscapine'; %maximize noscapine production rate
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
AA_Exchanges = {'r_1873' ... % L-alanine exchange
                'r_1879' ... % L-arginine exchange
                'r_1880' ... % L-asparagine exchange
                'r_1881' ... % L-aspartate exchange
                'r_1883' ... %	L-cysteine exchange
                'r_1889' ... %	L-glutamate exchange
                'r_1891' ... %	L-glutamine exchange
                'r_1810' ... %	L-glycine exchange
                'r_1893' ... %	L-histidine exchange
                'r_1897' ... %	L-isoleucine exchange
                'r_1899' ... %	L-leucine exchange
                'r_1900' ... %	L-lysine exchange
                'r_1902' ... %	L-methionine exchange
                'r_1903' ... %	L-phenylalanine exchange
                'r_1904' ... %	L-proline exchange
                'r_1906' ... %	L-serine exchange
                'r_1911' ... %	L-threonine exchange
                'r_1912' ... %	L-tryptophan exchange
                'r_1913' ... %	L-tyrosine exchange
                'r_1914'};    %	L-valine exchange
AA_IDs =   {'Ala' 'Arg' 'Asn' 'Asp' 'Cys' 'Glu' 'Gln' 'Gly' 'His' 'Ile' ...
            'Leu' 'Lys' 'Met' 'Phe' 'Pro' 'Ser' 'Thr' 'Trp' 'Tyr' 'Val'};

mu = 0.1;
model = changeRxnBounds(model,'r_2111',mu,'b');
q_aa = -0.01;

model_ref = model;
fileName = writeLP(model_ref,mu,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator);
command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
system(command,'-echo');
[~,~,flux_ref] = readSoplexResult('Simulation.lp.out',model_ref);
q_fe_ref = flux_ref(ismember(model_ref.rxnNames,'iron(3+) exchange'),1);

q_fe_1 = q_fe_ref;
model_ref_1 = model;
model_ref_1.lb(ismember(model_ref_1.rxnNames,'iron(3+) exchange')) = q_fe_1;
fileName = writeLP(model_ref_1,mu,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator);
command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
system(command,'-echo');
[~,~,flux_ref_1] = readSoplexResult('Simulation.lp.out',model_ref_1);

fluxes_aa_1 = zeros(length(model.rxns),length(AA_Exchanges));
for i = 1:length(AA_Exchanges)
    model_tmp = model_ref_1;
    model_tmp = changeRxnBounds(model_tmp,AA_Exchanges{i},q_aa,'l');
    disp(['Unlimited iron uptake: ' AA_IDs{i}]);
    fileName = writeLP(model_tmp,mu,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator);
    command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
    system(command,'-echo');
    [sol_obj,sol_status,sol_full] = readSoplexResult('Simulation.lp.out',model_tmp);
    disp(['solution status: ' sol_status]);
    if strcmp(sol_status,'optimal')
        fluxes_aa_1(:,i) = sol_full;
    end
end

q_fe_2 = q_fe_ref * 0.5;
model_ref_2 = model;
model_ref_2.lb(ismember(model_ref_2.rxnNames,'iron(3+) exchange')) = q_fe_2;
fileName = writeLP(model_ref_2,mu,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator);
command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
system(command,'-echo');
[~,~,flux_ref_2] = readSoplexResult('Simulation.lp.out',model_ref_2);

fluxes_aa_2 = zeros(length(model.rxns),length(AA_Exchanges));
for i = 1:length(AA_Exchanges)
    model_tmp = model_ref_2;
    model_tmp = changeRxnBounds(model_tmp,AA_Exchanges{i},q_aa,'l');
    disp(['50% iron uptake: ' AA_IDs{i}]);
    fileName = writeLP(model_tmp,mu,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator);
    command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
    system(command,'-echo');
    [sol_obj,sol_status,sol_full] = readSoplexResult('Simulation.lp.out',model_tmp);
    disp(['solution status: ' sol_status]);
    if strcmp(sol_status,'optimal')
        fluxes_aa_2(:,i) = sol_full;
    end
end

sN3_res.lables = ['Ref' AA_IDs];
sN3_res.fluxes_UnlimitedFe = [flux_ref_1 fluxes_aa_1];
sN3_res.fluxes_LimitedFe = [flux_ref_2 fluxes_aa_2];


save('sN3_res.mat','sN3_res');

clear;

toc;

