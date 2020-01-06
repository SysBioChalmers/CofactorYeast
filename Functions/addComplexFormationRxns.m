%% addComplexFormationRxns 
function model = addComplexFormationRxns(model,Protein_stoichiometry)

% find reactions with GPR
idx = ~ismember(model.grRules,'');
metrxnid_list = model.rxns(idx);
gpr_list = model.grRules(idx);

for i = 1:length(metrxnid_list)
    disp(['Adding complex formation:' num2str(i) '/' num2str(length(metrxnid_list))]);
    
    metrxnname = metrxnid_list(i);
    cmplxid = strcat(cell2mat(metrxnname),'_enzyme');
    rxnid = strcat(cell2mat(metrxnname),'_enzyme_formation');
    
    gpr = gpr_list(i);
    gpr_tmp = split(gpr,' and ');
    gpr_tmp_tmp = cellfun(@(x) strrep(x,'-','_'),gpr_tmp,'UniformOutput',false);
    subs_id = cellfun(@(x) strcat(x,'_cofactorbound'),gpr_tmp_tmp,'UniformOutput',false);
    
    % check if protein stoichiometry is determined
    idx_tmp = ismember(gpr_tmp,Protein_stoichiometry.protein);
    determined = gpr_tmp(idx_tmp);
    determined_subs_id = subs_id(idx_tmp);
    determined_stoich = Protein_stoichiometry.stoichiometry(ismember(Protein_stoichiometry.protein,determined));
    undetermined_subs_id = subs_id(~idx_tmp);
    undetermined_stoich = ones(length(undetermined_subs_id),1); % assumed to be 1 if not determined
    
    metlist = [determined_subs_id' undetermined_subs_id' cmplxid];
    coeflist = [-1*determined_stoich' -1*undetermined_stoich' 1];
    model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);
end