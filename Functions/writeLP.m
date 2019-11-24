%% writeLP 
function fileName = writeLP(model,mu,f,osenseStr,rxnID,enzymedata,factor_k)
% f is the fraction (g/gCDW) of the modeled proteins.

fileName = sprintf('Simulation.lp');
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
	if kcat < quantile(kcat_list,0.1,1)
        kcat = quantile(kcat_list,0.1,1);
	end
    kcat = kcat*factor_k;

    %find enzyme formation reaction id
    id_syn = strcat(enzyme,'_formation');
    idx_syn = find(strcmp(model.rxns,id_syn));
    
    %calculate the coefficient of enzyme formation reaction rate
    coef = kcat/mu;
    
    %find enzyme used reaction id (metabolic reaction)
    idx_rxn = find(ismember(model.rxns,strrep(enzyme,'_enzyme','')));

    fprintf(fptr,'CM%d: X%d - %.15f X%d <= 0\n',...
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
eq = sprintf('%s + %.15f X%d%c',eq,1,idx,sep);

fprintf(fptr,'Ctotprot: %s = %.15f\n',eq,mu*f);

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
