%% expandModel 
function [model,enzymedata] = expandModel(model,enzymedata)
% Assuming that kcat would be non-zero if no cofactors will be binding on
% the enzyme, but just a lower value than the maximum, the function here
% is to add new reactions for enzymes with cofactors to account for lower
% kcats cases.
disp('Expanding model...');

% add new enzyme formation and cofactor binding reactions for the enzymes
% that has cofactors, in which no cofactors will be added.
rxnidx_enzymeformation = find(contains(model.rxns,'_enzyme_formation'));
for i = 1:length(rxnidx_enzymeformation)
    disp([num2str(i) '/' num2str(length(rxnidx_enzymeformation))]);
    idx_tmp = rxnidx_enzymeformation(i);
    rxnid = model.rxns{idx_tmp};
    if ~contains(rxnid,'withoutcofactor') % check if the enzyme is the type without cofactor
        rxnid_check = strrep(rxnid,'_enzyme_formation','_withoutcofactor_enzyme_formation');
        if ~ismember(rxnid_check,model.rxns) % check if the enzyme already has the type without cofactor
            ezmid = strrep(rxnid,'_formation','');
            cf_count = sum(enzymedata.cofactor_copy(ismember(enzymedata.enzyme,ezmid),:));
            if cf_count > 0 % the enzyme has at least one cofactor
                % search for the cofactor binding reactions for its subunits
                subunit_id_list = model.mets(model.S(:,idx_tmp) < 0);
                subunit_s_list = full(model.S(ismember(model.mets,subunit_id_list),idx_tmp));
                for j = 1:length(subunit_id_list)
                    subunit_id = subunit_id_list(j);
                    cfbinding_rxn = strcat('r_',subunit_id);
                    if length(find(model.S(:,ismember(model.rxns,cfbinding_rxn)))) > 2
                        subs_tmp = strrep(subunit_id,'_cofactorbound','_translated');
                        prod_tmp = strcat(subunit_id,'_withoutcofactor');
                        rxnid_tmp = strcat('r_',prod_tmp);
                        model = addReaction(model,cell2mat(rxnid_tmp),'metaboliteList',[subs_tmp,prod_tmp],'stoichCoeffList',[-1,1],'reversible',false);
                        subunit_id_list(j) = prod_tmp;
                    end
                end
                % add new enzyme formation reaction
                rxn_ezmfmt_id = strrep(rxnid,'_enzyme_formation','_withoutcofactor_enzyme_formation');
                new_ezm_id = {strrep(ezmid,'_enzyme','_withoutcofactor_enzyme')};
                model = addReaction(model,rxn_ezmfmt_id,'metaboliteList',[subunit_id_list;new_ezm_id]','stoichCoeffList',[subunit_s_list;1]','reversible',false);
                % add new dilution reaction
                rxn_dil_id = strcat(new_ezm_id,'_dilution');
                model = addReaction(model,cell2mat(rxn_dil_id),'metaboliteList',new_ezm_id,'stoichCoeffList',-1,'reversible',false);
                % add new metabolic reaction as iso-reaction
                old_m_rxn_id = strrep(ezmid,'_enzyme','');
                m_idx = ismember(model.rxns,old_m_rxn_id);
                m_mets_list = model.mets(model.S(:,m_idx) ~= 0);
                m_mets_s_list = full(model.S(ismember(model.mets,m_mets_list),m_idx));
                m_gpr = model.grRules{m_idx};
                new_m_rxn_id = strcat(old_m_rxn_id,'_withoutcofactor');
                model = addReaction(model,new_m_rxn_id,'metaboliteList',m_mets_list',...
                        'stoichCoeffList',m_mets_s_list','geneRule',m_gpr,...
                        'lowerBound',model.lb(m_idx),'upperBound',model.ub(m_idx));
                % update enzyme data
                idx_ezm = ismember(enzymedata.enzyme,ezmid);
                enzymedata.enzyme = [enzymedata.enzyme;new_ezm_id];
                enzymedata.substrate = [enzymedata.substrate;enzymedata.substrate(idx_ezm,:)];
                enzymedata.subunit = [enzymedata.subunit;enzymedata.subunit(idx_ezm,:)];
                enzymedata.subunit_stoichiometry = [enzymedata.subunit_stoichiometry;enzymedata.subunit_stoichiometry(idx_ezm,:)];
                enzymedata.subunit_ec = [enzymedata.subunit_ec;enzymedata.subunit_ec(idx_ezm,:)];
                enzymedata.subunit_kcat = [enzymedata.subunit_kcat;enzymedata.subunit_kcat(idx_ezm,:)];
                enzymedata.subunit_kcat_conf = [enzymedata.subunit_kcat_conf;enzymedata.subunit_kcat_conf(idx_ezm,:)];
                enzymedata.kcat = [enzymedata.kcat;enzymedata.kcat(idx_ezm,:)];
                enzymedata.kcat_conf = [enzymedata.kcat_conf;enzymedata.kcat_conf(idx_ezm,:)];
                enzymedata.enzyme_MW = [enzymedata.enzyme_MW;enzymedata.enzyme_MW(idx_ezm,:)];
                enzymedata.cofactor_type = [enzymedata.cofactor_type;repelem({''},20)];
                enzymedata.cofactor_copy = [enzymedata.cofactor_copy;zeros(1,20)];
            end
        end
    end
end























