%% addCofactorRxns 
function model = addCofactorRxns(model)
% There are several assumptions:
% 1. 
% 2. 
% 3. no machineries included.

%% Import Cofactor Information
% [~,raw,~] = xlsread('xxx.xlsx');
bound_prot_list = cell(0,1);


%% Add reactions
% Collect all translated proteins
translated_prot_list = model.mets(contains(model.mets,'_translated'));

for i = 1:length(translated_prot_list)
    disp(['Adding cofactor binding:' num2str(i) '/' num2str(length(translated_prot_list))]);
    translated_prot_id = translated_prot_list{i};
    normal_id = strrep(translated_prot_id,'_translated','');
    normal_id = strrep(normal_id,'_','-');
    if ismember(normal_id,bound_prot_list)
        %%% add cofactor binding reactions
    else
        subsid = translated_prot_id;
        prodid = strrep(translated_prot_id,'_translated','_cofactorbound');
        rxnid = strcat('r_',prodid);
        metlist = [{subsid} {prodid}];
        coeflist = [-1 1];
        model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);
    end
end