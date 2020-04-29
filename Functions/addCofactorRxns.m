%% addCofactorRxns 
function model = addCofactorRxns(model)
% There are several assumptions:
% 1. 
% 2. 
% 3. no machineries included.

%% Import Cofactor Information
[~,raw,~] = xlsread('cofactorList.xlsx');
raw1 = raw(2:end,1);
raw2 = raw(2:end,2);
cofactorid = raw1(~ismember(raw2,''));
cofactormetid = raw2(~ismember(raw2,''));
load('CofactorDataset.mat');
idx = ismember(CofactorDataset.cofactor,cofactorid);
new_cofactor = CofactorDataset.cofactor(idx);
new_copy = CofactorDataset.copy(idx);
new_protein = CofactorDataset.protein(idx);

%% Add reactions
% Collect all translated proteins
translated_prot_list = model.mets(contains(model.mets,'_translated'));

for i = 1:length(translated_prot_list)
    disp(['Adding cofactor binding:' num2str(i) '/' num2str(length(translated_prot_list))]);
    translated_prot_id = translated_prot_list{i};
    normal_id = strrep(translated_prot_id,'_translated','');
    normal_id = strrep(normal_id,'_','-');
    if ismember(normal_id,new_protein)
        %%% add cofactor binding reactions
        idx_tmp = ismember(new_protein,normal_id);
        type = new_cofactor(idx_tmp);
        copy = new_copy(idx_tmp);
        subsid = translated_prot_id;
        prodid = strrep(translated_prot_id,'_translated','_cofactorbound');
        rxnid = strcat('r_',prodid);
        [~,b] = ismember(type,cofactorid);
        cf_metName = cofactormetid(b);
        cf_metid = model.mets(ismember(model.metNames,cf_metName));
        metlist = [{subsid} cf_metid' {prodid}];
        coeflist = [-1 -1*copy' 1];
        model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);
    else
        subsid = translated_prot_id;
        prodid = strrep(translated_prot_id,'_translated','_cofactorbound');
        rxnid = strcat('r_',prodid);
        metlist = [{subsid} {prodid}];
        coeflist = [-1 1];
        model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);
    end
end