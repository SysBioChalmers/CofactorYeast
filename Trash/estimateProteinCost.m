load('CofactorYeast.mat');
load('enzymedata.mat');

% block some reactions
model = blockRxns(model);

%% Produce co2
model_tmp = model;
model_tmp = changeObjective(model_tmp,'r_4046'); %
model_tmp = changeRxnBounds(model_tmp,'r_4046',0,'l');
model_tmp = changeRxnBounds(model_tmp,'r_4046',1000,'u');
% model_tmp = changeRxnBounds(model_tmp,'r_1672',6,'b');
sol = optimizeCbModel(model_tmp,'max','one');
fluxes = sol.x;
tot_protein_cost_co2 = calculateProteinCost(model_tmp,fluxes,enzymedata)/abs(fluxes(ismember(model_tmp.rxns,'r_1714')));
Yatp_co2 = fluxes(ismember(model_tmp.rxns,'r_4046'))/abs(fluxes(ismember(model_tmp.rxns,'r_1714')));
efficiency_co2 = fluxes(ismember(model_tmp.rxns,'r_4046'))/tot_protein_cost_co2;

%% Produce ethanol
model_tmp = model;
model_tmp = changeObjective(model_tmp,'r_4046'); %
model_tmp = changeRxnBounds(model_tmp,'r_4046',0,'l');
model_tmp = changeRxnBounds(model_tmp,'r_4046',1000,'u');
model_tmp = changeRxnBounds(model_tmp,'r_1761',2,'b');
sol = optimizeCbModel(model_tmp,'max','one');
fluxes = sol.x;
tot_protein_cost_ethanol = calculateProteinCost(model_tmp,fluxes,enzymedata)/abs(fluxes(ismember(model_tmp.rxns,'r_1714')));
Yatp_ethanol = fluxes(ismember(model_tmp.rxns,'r_4046'))/abs(fluxes(ismember(model_tmp.rxns,'r_1714')));
efficiency_ethanol = fluxes(ismember(model_tmp.rxns,'r_4046'))/tot_protein_cost_ethanol;

%% Produce acetaldehyde
model_tmp = model;
model_tmp = changeObjective(model_tmp,'r_4046'); %
model_tmp = changeRxnBounds(model_tmp,'r_4046',0,'l');
model_tmp = changeRxnBounds(model_tmp,'r_4046',1000,'u');
model_tmp = changeRxnBounds(model_tmp,'r_1631',2,'b');
sol = optimizeCbModel(model_tmp,'max','one');
fluxes = sol.x;
tot_protein_cost_acetaldehyde = calculateProteinCost(model_tmp,fluxes,enzymedata)/abs(fluxes(ismember(model_tmp.rxns,'r_1714')));
Yatp_acetaldehyde = fluxes(ismember(model_tmp.rxns,'r_4046'))/abs(fluxes(ismember(model_tmp.rxns,'r_1714')));
efficiency_acetaldehyde = fluxes(ismember(model_tmp.rxns,'r_4046'))/tot_protein_cost_acetaldehyde;

%% Produce acetate
model_tmp = model;
model_tmp = changeObjective(model_tmp,'r_4046'); %
model_tmp = changeRxnBounds(model_tmp,'r_4046',0,'l');
model_tmp = changeRxnBounds(model_tmp,'r_4046',1000,'u');
model_tmp = changeRxnBounds(model_tmp,'r_1634',2,'b');
sol = optimizeCbModel(model_tmp,'max','one');
fluxes = sol.x;
tot_protein_cost_acetate = calculateProteinCost(model_tmp,fluxes,enzymedata)/abs(fluxes(ismember(model_tmp.rxns,'r_1714')));
Yatp_acetate = fluxes(ismember(model_tmp.rxns,'r_4046'))/abs(fluxes(ismember(model_tmp.rxns,'r_1714')));
efficiency_acetate = fluxes(ismember(model_tmp.rxns,'r_4046'))/tot_protein_cost_acetate;


