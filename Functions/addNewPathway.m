%% addCofactorRxns
function [model,enzymedata] = addNewPathway(model,enzymedata,file_name)

% load data
[~,rxnlst,~] = xlsread(file_name,'Metabolic_rxns');
rxnlst = rxnlst(2:end,:);

[~,gprlst,~] = xlsread(file_name,'Gpr_info');
[gpr_cfclst,~,~] = xlsread(file_name,'Gpr_info');
gpr_gprlst = strtrim(gprlst(2:end,1));
gpr_seqlst = strtrim(gprlst(2:end,3));
gpr_cftlst = strtrim(gprlst(2:end,5));

[~,enzymelst,~] = xlsread(file_name,'Enzyme_info');
[ezm_gprslst,~,~] = xlsread(file_name,'Enzyme_info');
ezm_idlst = strtrim(enzymelst(2:end,1));
ezm_gprlst = strtrim(enzymelst(2:end,2));

[~,kcatlst,~] = xlsread(file_name,'kcat_info');
[kcats,~,~] = xlsread(file_name,'kcat_info');
kcat_rxnlst = strtrim(kcatlst(2:end,1));
kcats = kcats * 3600; % unit: /h


%% add metabolic rxns
for i = 1:length(rxnlst(:,1))
    disp(['Adding metabolic reactions:' num2str(i) '/' num2str(length(rxnlst(:,1)))]);
    rxnformula = strrep(rxnlst{i,3},' [','[');
    [metaboliteList, stoichCoeffList, revFlag] = parseRxnFormula(rxnformula);
    metaboliteList = strrep(metaboliteList,'[',' [');
    metaboliteList = strrep(metaboliteList,'&',' ');
    comps = split(metaboliteList', ' [');
    if length(comps(1,:)) == 1
        comps = comps(2);
    else
        comps = comps(:,2);
    end
    comps = strrep(comps,']','');
    CONValldata = cat(2,model.compNames,model.comps);
    [~,b] = ismember(comps,CONValldata(:,1));
    comps = CONValldata(b,2);
    
    %mapping mets to model.metnames, get s_ index for new mets
    for j = 1:length(metaboliteList)
        [~,metindex] = ismember(metaboliteList(j),model.metNames);
        if metindex ~= 0
            mets(j) = model.mets(metindex);
        elseif metindex == 0
            newID = sum(startsWith(model.mets,'new_m_')) + 1;
            newID = num2str(newID);
            mets(j) = strcat('new_m_',newID,'[',comps(j),']');
            model = addMetabolite(model,char(mets(j)), ...
                'metName',metaboliteList{j});
        end
    end
    
    model = addReaction(model,...
        rxnlst{i,1},...
        'reactionName', rxnlst{i,2},...
        'metaboliteList',mets,...
        'stoichCoeffList',stoichCoeffList,...
        'reversible',revFlag,...
        'geneRule',rxnlst{i,4},...
        'checkDuplicate',1);
    clear mets;
end

%% add translation and cofactor binding rxns

[~,raw,~] = xlsread('aa_id.xlsx');
aa_list = struct();
aa_list.aa = raw(2:end,1);
aa_list.subs = raw(2:end,3);
aa_list.prod = raw(2:end,5);

for i = 1:length(gpr_gprlst)
    disp(['Adding translation and cofactor binding reactions:' num2str(i) '/' num2str(length(gpr_gprlst))]);
    gprid = gpr_gprlst(i);
    seq = gpr_seqlst(i);
    tot = countAA(seq,aa_list);
    protid = cell2mat(gprid);
    protid = strcat(protid,'_translated');
    rxnid = strcat('r_',protid);
    metlist = [tot.subs' tot.prod' protid];
    coeflist = [-1*tot.num' tot.num' 1];
    model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);
    
    subsid = protid;
    prodid = strrep(protid,'_translated','_cofactorbound');
    rxnid = strcat('r_',prodid);
    if ismember(gpr_cftlst(i),'-')
        metlist = [{subsid} {prodid}];
        coeflist = [-1 1];
    else
        if ismember(gpr_cftlst(i),'FE_II')
            cf_metName = 'iron(2+) [mitochondrion]';
            [~,c] = ismember(cf_metName,model.metNames);
            cf_metid = model.mets(c);
            cfc_tmp = gpr_cfclst(i);
            metlist = [{subsid} cf_metid {prodid}];
            coeflist = [-1 -1*cfc_tmp 1];
        elseif ismember(gpr_cftlst(i),'HEME_A')
            cf_metName = 'heme a [mitochondrion]';
            [~,c] = ismember(cf_metName,model.metNames);
            cf_metid = model.mets(c);
            cfc_tmp = gpr_cfclst(i);
            metlist = [{subsid} cf_metid {prodid}];
            coeflist = [-1 -1*cfc_tmp 1];
        end
    end
    model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);
end

%% add complex formation and dilution rxns
unq_ezm_list = unique(ezm_idlst);

for i = 1:length(unq_ezm_list)
    disp(['Adding complex formation and dilution reactions:' num2str(i) '/' num2str(length(unq_ezm_list))]);
    ezm_tmp = unq_ezm_list(i);
    ezm_id_tmp = strcat(ezm_tmp,'_enzyme');
    rxn_fmt_id = strcat(ezm_tmp,'_enzyme_formation');
    rxn_dil_id = strcat(ezm_tmp,'_enzyme_dilution');
    
    subunit_tmp = ezm_gprlst(ismember(ezm_idlst,ezm_tmp));
    subunit_s_tmp = ezm_gprslst(ismember(ezm_idlst,ezm_tmp));
    metlist = [subunit_tmp' ezm_id_tmp];
    coeflist = [-1*subunit_s_tmp' 1];
    model = addReaction(model,cell2mat(rxn_fmt_id),'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);
    
    model = addReaction(model,cell2mat(rxn_dil_id),'metaboliteList',ezm_id_tmp,'stoichCoeffList',-1,'reversible',false);
end

%% add enzyme information
for i = 1:length(kcat_rxnlst)
    disp(['Adding enzyme information:' num2str(i) '/' num2str(length(kcat_rxnlst))]);
    ezm_tmp = kcat_rxnlst(i);
    ezm_id_tmp = strcat(ezm_tmp,'_enzyme');
    enzymedata.enzyme = [enzymedata.enzyme; ezm_id_tmp];
    enzymedata.substrate = [enzymedata.substrate; {'skip'}];
    
    subunit_tmp = ezm_gprlst(ismember(ezm_idlst,ezm_tmp));
    subunit_s_tmp = ezm_gprslst(ismember(ezm_idlst,ezm_tmp));
    enzymedata.subunit = [enzymedata.subunit; [subunit_tmp' repelem({''},20-length(subunit_tmp))]];
    enzymedata.subunit_stoichiometry = [enzymedata.subunit_stoichiometry; [subunit_s_tmp' zeros(1,20-length(subunit_s_tmp))]];
    
    enzymedata.subunit_ec = [enzymedata.subunit_ec; repelem({'skip'},20)];
    enzymedata.subunit_kcat = [enzymedata.subunit_kcat; nan(1,20)];
    enzymedata.subunit_kcat_conf = [enzymedata.subunit_kcat_conf; nan(1,20)];
    
    enzymedata.kcat = [enzymedata.kcat; kcats(i)];
    enzymedata.kcat_conf = [enzymedata.kcat_conf; nan];
    
    MW_sum = 0;
    for j = 1:length(subunit_tmp)
        subunit_id_tmp = subunit_tmp(j);
        subunit_seq_tmp = gpr_seqlst(ismember(gpr_gprlst,subunit_id_tmp));
        subunit_MW_tmp = getProteinMW(cell2mat(subunit_seq_tmp));
        MW_sum = MW_sum + subunit_MW_tmp * subunit_s_tmp(j);
    end
    enzymedata.enzyme_MW = [enzymedata.enzyme_MW; MW_sum];
    
    type = cell(1,0);
    copy = zeros(1,0);
    for j = 1:length(subunit_tmp)
        subunit_id_tmp = subunit_tmp(j);
        rxn_idx_tmp = ismember(model.rxns,strcat('r_',subunit_id_tmp,'_cofactorbound'));
        met_idx_tmp = model.S(:,rxn_idx_tmp) ~= 0;
        metlist = model.mets(met_idx_tmp);
        slist = full(model.S(met_idx_tmp,rxn_idx_tmp));
        typelist = metlist(~contains(metlist,subunit_id_tmp));
        copylist = -1*slist(~contains(metlist,subunit_id_tmp));
        type = [type typelist'];
        copy = [copy copylist' * subunit_s_tmp(j)];
    end
    uniq_type = unique(type);
    tmp_type = cell(1,20);
    tmp_copy = zeros(1,20);
    if ~isempty(uniq_type)
        for k = 1:length(uniq_type)
            tmp_type(1,k) = uniq_type(k);
            tmp_copy(1,k) = sum(copy(ismember(type,uniq_type(k))));
        end
    end
    enzymedata.cofactor_type = [enzymedata.cofactor_type; tmp_type];
    enzymedata.cofactor_copy = [enzymedata.cofactor_copy; tmp_copy];
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MW = getProteinMW(seq)
% (From GECKO)
aa_codes = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N', ...
            'O','P','Q','R','S','T','U','V','W','X','Y','Z'};
aa_MWs   = [71.08 114.60 103.14 115.09 129.11 147.17 57.05 137.14 ...
            113.16 113.16 128.17 113.16 131.20 114.10 255.31 97.12 ...
            128.13 156.19 87.08 101.10 150.04 99.13 186.21 126.50 ...
            163.17 128.62];
MW = 18;
for i = 1:length(aa_codes)
    count = length(strfind(seq,aa_codes{i}));
    MW = MW + count*aa_MWs(i);
end
end