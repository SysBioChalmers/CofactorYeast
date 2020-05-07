%% simulationGeneKnockout
function simulationGeneKnockout(a,b)
Startme
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
gene_list = model.genes;
% load('sGE_tf.mat');
% gene_list = model.genes(~sGE_tf);

fluxes = zeros(length(model.rxns),b-a+1);
genes = cell(1,b-a+1);

for i = a:b
    geneid_tmp = gene_list{i};
    geneid_tmp = strrep(geneid_tmp,'-','_');
    rxnid_tmp = strcat('r_',geneid_tmp,'_translated');
    model_tmp = changeRxnBounds(model,rxnid_tmp,0,'b');% block translation
    [~,flux_tmp] = searchMaxgrowthSpecial(model_tmp,geneid_tmp,f,osenseStr,rxnID,enzymedata,1e-6);
    fluxes(:,i-a+1) = flux_tmp;
    genes(1,i-a+1) = gene_list(i);
end

filename_1 = strcat('sGK_fluxes_',num2str(a),'.mat');
filename_2 = strcat('sGK_genes_',num2str(a),'.mat');
cd Results/;
save(filename_1,'fluxes');
save(filename_2,'genes');
cd ../;
end
