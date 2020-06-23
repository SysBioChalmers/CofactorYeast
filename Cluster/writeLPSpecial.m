%% writeLP 
function fileName = writeLPSpecial(model,label_tmp,mu,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator,factor_k)
% f is the fraction (g/gCDW) of the modeled proteins.
% f_mito is the fraction (g/gCDW) of the mitochondrial proteins.
if exist('factor_k', 'var')
    if isempty(factor_k)
        factor_k = 1;
    end
else
    factor_k = 1;
end

fileNametmp = strcat('Simulation_',label_tmp,'.lp');
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

% newly added pathway reactions
newidx = contains(enzyme_list,'new_r_');
expandidx = contains(enzyme_list,'withoutcofactor');
lowkcat = quantile(kcat_list(~(newidx|expandidx)),0.05,1);

for i = 1:length(enzyme_list)
    enzyme = enzyme_list{i};
  	kcat = kcat_list(i);
    
    % Change kcats extremely low for original enzymes
    if ismember(enzyme,enzyme_list(~newidx))
        if kcat < lowkcat
            kcat = lowkcat;
        end
    end
	if contains(enzyme,'withoutcofactor')
        kcat = kcat*factor_k_withoutcofator;
	end
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
% 4) Constraint on mitochondrial proteins.
[~, mitoRxns] = collectCompartment(model,{'m','mm'});
mito_dil_rxns = strcat(mitoRxns,'_enzyme_dilution');

for i = 1:length(mito_dil_rxns)
    rxn_id = mito_dil_rxns{i};
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

fprintf(fptr,'Cmito: %s <= %.15f\n',eq,mu*f_mito);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5) Dilution of cofactors on the dummy complex

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