%% addTranslationRxns 
function model = addTranslationRxns(model,ProteinSequence)
% There are several assumptions:
% 1. all proteins are translated in cytoplasm;
% 2. no energy cost.
% 3. no machineries included, e.g., ribosomes.

%% Import AA IDs
[~,raw,~] = xlsread('aa_id.xlsx');
aa_list = struct();
aa_list.aa = raw(2:end,1);
aa_list.subs = raw(2:end,3);
aa_list.prod = raw(2:end,5);

%% Add reactions
for i = 1:length(model.genes)
    disp(['Adding translation:' num2str(i) '/' num2str(length(model.genes))]);
    geneid = model.genes(i);
    seq = ProteinSequence.seq(ismember(ProteinSequence.id,geneid));
    sum = countAA(seq,aa_list);
    protid = cell2mat(geneid);
    protid = strrep(protid,'-','_');
    protid = strcat(protid,'_translated');
    rxnid = strcat('r_',protid);
    metlist = [sum.subs' sum.prod' protid];
    coeflist = [-1*sum.num' sum.num' 1];
    model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);
end