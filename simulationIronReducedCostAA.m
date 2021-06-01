%% simulationIronReducedCostAA
% Add each amino acid to see how growth changes

% Timing: ~ 50000 s
tic;
load('CofactorYeast.mat');
load('enzymedata.mat');

soplexpath = '/Users/cheyu/build/bin/soplex'; % change this to the soplex path on your PC


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

%% Solve LPs
% reference
[~,flux_ref100] = searchMaxgrowth(model,f,f_mito,osenseStr,rxnID,enzymedata,0,1e-6,soplexpath);
q_fe_ref100 = flux_ref100(strcmp(model.rxns,'r_1861'),1);

%% Solve LPs

factor_k_withoutcofator = 0.5;
lower_fe = 0.5;
q_fe_new = q_fe_ref100*lower_fe;

model_fe50 = changeRxnBounds(model,'r_1861',q_fe_new,'b');
[~,flux_ref50] = searchMaxgrowth(model_fe50,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator,1e-6,soplexpath);


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
q_aa = -0.01;

fluxes = zeros(length(model.rxns),length(AA_Exchanges));
for i = 1:length(AA_Exchanges)
    model_tmp = model_fe50;
    model_tmp = changeRxnBounds(model_tmp,AA_Exchanges{i},q_aa,'l');
    disp(['Adding ' AA_IDs{i}]);
    [~,flux_tmp] = searchMaxgrowth(model_tmp,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator,1e-6,soplexpath);
    fluxes(:,i) = flux_tmp;
end

sIRCAA_res.lables = ['Ref' AA_IDs];
sIRCAA_res.fluxes = [flux_ref50 fluxes];

cd Results/;
save('sIRCAA_res.mat','sIRCAA_res');
cd ../;
clear;

toc;

