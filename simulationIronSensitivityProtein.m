%% simulationIronSensitivityProtein
% Protein level
% Assume that an iron-containing enzyme does not need iron for activity.
% 
% Timing: ~ 55000 s
tic;
load('CofactorYeast.mat');
load('enzymedata.mat');

% collect iron-containing proteins in the model
iron_metNames = {'iron(3+) [mitochondrion]';
                 'iron(2+) [mitochondrion]';
                 'heme a [mitochondrion]';
                 'heme c [mitochondrion]';
                 'ferroheme b [mitochondrion]';
                 'Fe2S2 cluster [mitochondrion]';
                 'Fe4S4 cluster [mitochondrion]';
                 'Fe3S4 cluster [mitochondrion]'};
iron_mets = model.mets(ismember(model.metNames,iron_metNames));

iron_containing_proteins_tmp = cell(1,0);
for i = 1:length(model.genes)
%     disp([num2str(i) '/' num2str(length(model.genes))]);
    cfb_rxn = strcat('r_',strrep(model.genes(i),'-','_'),'_cofactorbound');
    if ismember(cfb_rxn,model.rxns)
        rxnidx = ismember(model.rxns,cfb_rxn);
        cflist = model.mets(model.S(:,rxnidx) < 0);
        if any(contains(cflist,iron_mets))
            iron_containing_proteins_tmp = [iron_containing_proteins_tmp;model.genes(i)];
        end
    end
end

% remove iron-containing enzymes without fluxes according to previous
% simulations
load('sI_res.mat');
pre_fluxes = [sI_res.flux_ref sI_res.fluxes];
nonzero_rxns = model.rxns(any(pre_fluxes,2));
iron_containing_proteins = cell(1,0);
for i = 1:length(iron_containing_proteins_tmp)
    protein_tmp = iron_containing_proteins_tmp(i);
    protein_tmp_tmp = [strcat('r_',strrep(protein_tmp,'-','_'),'_cofactorbound');
                       strcat('r_',strrep(protein_tmp,'-','_'),'_cofactorbound_withoutcofactor')];
    if any(ismember(protein_tmp_tmp,nonzero_rxns))
        iron_containing_proteins = [iron_containing_proteins;protein_tmp];
    end
end

clear cofactor_tmp protein_tmp protein_tmp_tmp i iron_containing_proteins_tmp;
clear iron_metNames nonzero_rxns pre_fluxes sI_res rxnidx cflist cfb_rxn;

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
[~,flux_ref] = searchMaxgrowth(model,f,f_mito,osenseStr,rxnID,enzymedata,0,1e-6);
q_fe_ref = flux_ref(strcmp(model.rxns,'r_1861'),1);

%% Solve LPs

factor_k_withoutcofator = 0.5;
lower_fe = 0.5;

fluxes = zeros(length(model.rxns),length(iron_containing_proteins));
for i = 1:length(iron_containing_proteins)
    protein = iron_containing_proteins(i);
    disp([num2str(i) '/' num2str(length(iron_containing_proteins))]);
    model_tmp = model;
    model_tmp = changeRxnBounds(model_tmp,'r_1861',q_fe_ref*lower_fe,'b');
    % block without reaction
    model_tmp.ub(ismember(model_tmp.rxns,strcat('r_',strrep(protein,'-','_'),'_cofactorbound_withoutcofactor'))) = 0;
    model_tmp.lb(ismember(model_tmp.rxns,strcat('r_',strrep(protein,'-','_'),'_cofactorbound_withoutcofactor'))) = 0;
    % remove iron
    rxn_cfb = strcat('r_',strrep(protein,'-','_'),'_cofactorbound');
    idx_cfb = ismember(model_tmp.rxns,rxn_cfb);
    idx_fe = model_tmp.S(:,idx_cfb) < 0 & ismember(model_tmp.mets,iron_mets);
    model_tmp.S(idx_fe,idx_cfb) = 0;
    
    [~,flux_tmp] = searchMaxgrowth(model_tmp,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator,1e-6);
    fluxes(:,i) = flux_tmp;
end

sISP_res.proteins = iron_containing_proteins;
sISP_res.fluxes = fluxes;
sISP_res.flux_ref = flux_ref;

cd Results/;
save('sISP_res.mat','sISP_res');
cd ../;
clear;

toc;

