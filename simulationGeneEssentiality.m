%% simulationGeneEssentiality
% Timing: ~ 65000 s
tic;
load('CofactorYeast.mat');
load('enzymedata.mat');

%% Set model
% set medium
model = setMedia(model,3);% YNB+CSM-Ura
% set carbon source
model = changeRxnBounds(model,'r_1714',-1000,'l');% glucose
% set oxygen
model = changeRxnBounds(model,'r_1992',-1000,'l');
% block reactions
model = blockRxns(model);
% model = changeRxnBounds(model,'r_1631',0,'b');% acetaldehyde production

%% Set optimization
rxnID = 'dilute_dummy';
osenseStr = 'Maximize';

tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
clear tot_protein f_modeled_protein;

%% Solve LPs
% reference
disp('Simulating reference state... ');
[mu_ref,~] = searchMaxgrowth(model,f,osenseStr,rxnID,enzymedata,1e-6);
disp(['Reference growth rate: ' num2str(round(mu_ref,4)) ' /h']);

% filter out essential genes
mu_cutoff = mu_ref * 0.01; % cutoff of essential genes
model_cutoff = changeRxnBounds(model,'r_2111',mu_cutoff,'b');

sGE_tf = false(length(model.genes),1);

for i = 1:length(model.genes)
    geneid_tmp = model.genes{i};
    disp(['Gene essentiality analysis ' num2str(i) '/' num2str(length(model.genes)) ': ' num2str(geneid_tmp)]);
    geneid_tmp = strrep(geneid_tmp,'-','_');
    rxnid_tmp = strcat('r_',geneid_tmp,'_translated');
    model_tmp = changeRxnBounds(model_cutoff,rxnid_tmp,0,'b');% block translation
    fileName = writeLP(model_tmp,mu_cutoff,f,osenseStr,rxnID,enzymedata,1);
    command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -t300 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
    system(command,'-echo');
    [~,sol_status,~] = readSoplexResult('Simulation.lp.out',model_tmp);
    disp(['solution status: ' sol_status]);
    sGE_tf(i) = ~strcmp(sol_status,'optimal');
end

cd Results/;
save('sGE_tf.mat','sGE_tf');
cd ../;
clear;

toc;


