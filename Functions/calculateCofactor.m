%% calculateCofactor 
function enzymedata = calculateCofactor(enzymedata,model)

enzymedata.cofactor_type = cell(length(enzymedata.enzyme),20);
enzymedata.cofactor_copy = zeros(length(enzymedata.enzyme),20);

for i = 1:length(enzymedata.enzyme)
    subunit_list = enzymedata.subunit(i,:);
    subunit_list = subunit_list(~ismember(subunit_list,''));
    type = cell(1,0);
    copy = zeros(1,0);
    for j = 1:length(subunit_list) 
        subunit_tmp = subunit_list(j);
        rxn_idx_tmp = ismember(model.rxns,strcat('r_',subunit_tmp,'_cofactorbound'));
        met_idx_tmp = model.S(:,rxn_idx_tmp) ~= 0;
        metlist = model.mets(met_idx_tmp);
        slist = full(model.S(met_idx_tmp,rxn_idx_tmp));
        typelist = metlist(~contains(metlist,subunit_tmp));
        copylist = -1*slist(~contains(metlist,subunit_tmp));
        type = [type typelist'];
        copy = [copy copylist' * enzymedata.subunit_stoichiometry(i,j)];
    end
    uniq_type = unique(type);
    if ~isempty(uniq_type)
        for k = 1:length(uniq_type)
            enzymedata.cofactor_type(i,k) = uniq_type(k);
            enzymedata.cofactor_copy(i,k) = sum(copy(ismember(type,uniq_type(k))));
        end
    end
end


