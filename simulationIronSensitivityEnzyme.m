%% simulationIronSensitivityEnzyme
% Enzyme level
% Assume that an iron-containing enzyme does not need iron for activity.
% 
% Timing: ~ 74000 s
tic;
load('CofactorYeast.mat');
load('enzymedata.mat');
load('CofactorDataset.mat');

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

iron_containing_enzymes_tmp = cell(1,0);
for i = 1:length(enzymedata.enzyme)
%     disp([num2str(i) '/' num2str(length(enzymedata.enzyme))]);
    
    cofactor_tmp = enzymedata.cofactor_type(i,:);
    cofactor_tmp = cellfun(@(x) num2str(x),cofactor_tmp,'UniformOutput',false);
    cofactor_tmp = cofactor_tmp(~ismember(cofactor_tmp,''));
    
    if any(contains(cofactor_tmp,iron_mets)) && ~contains(enzymedata.enzyme(i),'withoutcofactor')
        iron_containing_enzymes_tmp = [iron_containing_enzymes_tmp;enzymedata.enzyme(i)];
    end
end

% remove iron-containing enzymes without fluxes according to previous
% simulations
load('sI_res.mat');
pre_fluxes = [sI_res.flux_ref sI_res.fluxes];
nonzero_rxns = model.rxns(any(pre_fluxes,2));
iron_containing_enzymes = cell(1,0);
for i = 1:length(iron_containing_enzymes_tmp)
    enzymes_tmp = iron_containing_enzymes_tmp(i);
    enzymes_tmp_tmp = [strrep(enzymes_tmp,'_enzyme','');strrep(enzymes_tmp,'_enzyme','_withoutcofactor')];
    if any(ismember(enzymes_tmp_tmp,nonzero_rxns))
        iron_containing_enzymes = [iron_containing_enzymes;enzymes_tmp];
    end
end

clear cofactor_tmp enzymes_tmp enzymes_tmp_tmp i iron_containing_enzymes_tmp;
clear iron_metNames nonzero_rxns pre_fluxes sI_res;

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

fluxes = zeros(length(model.rxns),length(iron_containing_enzymes));
for i = 1:length(iron_containing_enzymes)
    enzyme = iron_containing_enzymes(i);
    disp([num2str(i) '/' num2str(length(iron_containing_enzymes))]);
    model_tmp = model;
    model_tmp = changeRxnBounds(model_tmp,'r_1861',q_fe_ref*lower_fe,'b');
    % block without reaction
    model_tmp.ub(ismember(model_tmp.rxns,strrep(enzyme,'_enzyme','_withoutcofactor'))) = 0;
    model_tmp.lb(ismember(model_tmp.rxns,strrep(enzyme,'_enzyme','_withoutcofactor'))) = 0;
    % remove iron, i.e. add iron as products for the enzyme formation reaction
    enzyme_fmt = strcat(enzyme,'_formation');
    cft_type = enzymedata.cofactor_type(ismember(enzymedata.enzyme,enzyme),:);
    cft_copy = enzymedata.cofactor_copy(ismember(enzymedata.enzyme,enzyme),:);
    cft_type = cellfun(@(x) num2str(x),cft_type,'UniformOutput',false);
    cft_copy = cft_copy(cft_copy ~= 0);
    cft_type = cft_type(~ismember(cft_type,''));
    cft_copy = cft_copy(ismember(cft_type,iron_mets));
    cft_type = cft_type(ismember(cft_type,iron_mets));
    [~,b] = ismember(cft_type,model_tmp.mets);
    model_tmp.S(b,ismember(model_tmp.rxns,enzyme_fmt)) = cft_copy;
    
    [~,flux_tmp] = searchMaxgrowth(model_tmp,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator,1e-6);
    fluxes(:,i) = flux_tmp;
end

sISE_res.enzymes = iron_containing_enzymes;
sISE_res.fluxes = fluxes;
sISE_res.flux_ref = flux_ref;

cd Results/;
save('sISE_res.mat','sISE_res');
cd ../;
clear;

toc;

