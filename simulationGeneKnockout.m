%% simulationGeneKnockout
function simulationGeneKnockout(a,b)
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
    fluxes(:,i) = flux_tmp;
    genes(1,i) = gene_list(i);
end

filename_1 = strcat('sGK_fluxes_',num2str(a),'.mat');
filename_2 = strcat('sGK_genes_',num2str(a),'.mat');
cd Results/;
save(filename_1,'fluxes');
save(filename_2,'genes');
cd ../;
end

function [mu,fluxes] = searchMaxgrowthSpecial(model,geneid_tmp,f,osenseStr,rxnID,enzymedata,precision,factor_k)

if exist('factor_k', 'var')
    if isempty(factor_k)
        factor_k = 1;
    end
else
    factor_k = 1;
end

if exist('precision', 'var')
    if isempty(precision)
        precision = 0.001;
    end
else
    precision = 0.001;
end

mu_low = 0;
mu_high = 1;

while mu_high-mu_low > precision
    mu_mid = (mu_low+mu_high)/2;
    model_tmp = changeRxnBounds(model,'r_2111',mu_mid,'b');
    disp(['mu = ' num2str(mu_mid)]);
    fileName = writeLPSpecial(model_tmp,geneid_tmp,mu_mid,f,osenseStr,rxnID,enzymedata,factor_k);
    command = sprintf('/cephyr/users/feiranl/Hebbe/tools/build/bin/soplex -s0 -g5 -t300 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
    system(command,'-echo');
    fileName_tmp = strcat('Simulation_',geneid_tmp,'.lp.out');
    [~,sol_status,sol_full] = readSoplexResult(fileName_tmp,model_tmp);
    disp(['solution status: ' sol_status]);
    
    if strcmp(sol_status,'optimal')
        mu_low = mu_mid;
        flux_tmp = sol_full;
    else
        mu_high = mu_mid;
    end
end

if mu_low > precision
    mu = mu_low;
    fluxes = flux_tmp;
else
    mu = 0;
    fluxes = zeros(length(model.rxns),1);
end
end

%% writeLP 
function fileName = writeLPSpecial(model,geneid_tmp,mu,f,osenseStr,rxnID,enzymedata,factor_k)
% f is the fraction (g/gCDW) of the modeled proteins.

fileNametmp = strcat('Simulation_',geneid_tmp,'.lp')
fileName = sprintf(fileNametmp);
fptr = fopen(fileName,'w');

% Set objective function
osenseStr = strcat(osenseStr,'\n');
fprintf(fptr,osenseStr);
index_obj = find(ismember(model.rxns,rxnID));
fprintf(fptr,'obj: X%d\n',index_obj);
fprintf(fptr,'Subject To\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) SV = 0.
% (From Ibrahim)
for i = 1:numel(model.mets)
    j = find(full(model.S(i,:)));
    for m = 1:numel(j)
        s = full(model.S(i,j(m)));
        if mod(m,200) == 0
            sep = newline;
        else
            sep = '';
        end
        if m == 1
           eq = sprintf('%.15f X%d',s,j(m));
        else
           if s>0
               eq = sprintf('%s + %.15f X%d%c',eq,s,j(m),sep);
           else
               eq = sprintf('%s %.15f X%d%c',eq,s,j(m),sep);
           end
        end
    end
    fprintf(fptr,'C%d: %s = 0\n',i,eq);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Coupling metabolic reactions and enzymes.

enzyme_list = enzymedata.enzyme;
kcat_list = enzymedata.kcat;


for i = 1:length(enzyme_list)
    enzyme = enzyme_list{i};
  	kcat = kcat_list(i);
    
% Change kcats extremely low
	if kcat < quantile(kcat_list,0.05,1)
        kcat = quantile(kcat_list,0.05,1);
	end
% 	if kcat < 3600
%         kcat = 3600;
% 	end
    
    kcat = kcat*factor_k;

    %find enzyme formation reaction id
    id_syn = strcat(enzyme,'_formation');
    idx_syn = find(strcmp(model.rxns,id_syn));
    
    %calculate the coefficient of enzyme formation reaction rate
    coef = kcat/mu;
    
    %find enzyme used reaction id (metabolic reaction)
    idx_rxn = find(ismember(model.rxns,strrep(enzyme,'_enzyme','')));

%     fprintf(fptr,'CM%d: X%d - %.15f X%d <= 0\n',...
%                  i,idx_rxn,coef,idx_syn);
    fprintf(fptr,'CM%d: X%d - %.15f X%d = 0\n',...
                 i,idx_rxn,coef,idx_syn);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Constraint on total enzymes and dummy complex

% total enzymes
dil_rxns = model.rxns(contains(model.rxns,'_dilution'));
for i = 1:length(dil_rxns)
    rxn_id = dil_rxns{i};
    comp_name = strrep(rxn_id,'_dilution','');
    idx = find(strcmp(model.rxns,rxn_id));
    
    if mod(i,200) == 0
        sep = newline;
	else
        sep = '';
    end

	MW = enzymedata.enzyme_MW(contains(enzymedata.enzyme,comp_name));
	coeff = MW/1000;
	if i == 1
        eq = sprintf('%.15f X%d',coeff,idx);
	else
        eq = sprintf('%s + %.15f X%d%c',eq,coeff,idx,sep);
	end
end

% dummy complex
idx = find(strcmp(model.rxns,'dilute_dummy'));
eq = sprintf('%s + %.15f X%d%c',eq,460/1000,idx,sep); %460 is MW of dummy complex (g/mol)

fprintf(fptr,'Ctotprot: %s = %.15f\n',eq,mu*f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4) Dilution of cofactors on the dummy complex

idx_dmycplx = find(strcmp(model.rxns,'dilute_dummy'));
idx_dmycofactor = find(strcmp(model.rxns,'dilute_dummy_cofactor'));
coef = 1000*f/460;
fprintf(fptr,'Cdummycofactor: X%d - %.15f X%d = 0\n',...
    idx_dmycplx,coef,idx_dmycofactor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set lower and upper bounds.

fprintf(fptr,'Bounds\n');

for i = 1:length(model.rxns)
	if model.ub(i) >= 100
        fprintf(fptr,'%f <= X%d <= +infinity\n',model.lb(i),i);
    else
        fprintf(fptr,'%f <= X%d <= %f\n',model.lb(i),i,model.ub(i));
	end
end

fprintf(fptr,'End\n');
fclose(fptr);
end


