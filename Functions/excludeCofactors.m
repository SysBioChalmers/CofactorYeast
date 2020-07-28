%% Exclude cofactors
function excludeCofactors
% Some enzymes are reported to contain cofactor just because they are
% catalysts of the synthesise of the cofactor. Therefore, here we will find
% out all the enzymes whose reactions contain the corresponding cofactor as
% substrates or products.
% These enzymes should be excluded for the cofactor.
load('Yeast8.mat');
[~,raw,~] = xlsread('cofactorList.xlsx');
raw1 = raw(2:end,1);
raw2 = raw(2:end,2);
cofactorid = raw1(~ismember(raw2,''));
cofactormetid = raw2(~ismember(raw2,''));

model = buildRxnGeneMat(org_model);
excludedata.cofactorid = cofactorid;
excludedata.geneid = cell(length(cofactorid),1);
for i = 1:length(cofactorid)
    metid = cofactormetid{i};
    metid = strsplit(metid,' [');
    metid = metid(1);
    metlist = strcat(metid,' [',model.compNames,']');
    rxnidx = any(model.S(ismember(model.metNames,metlist),:),1);
    genelist = model.genes(any(model.rxnGeneMat(rxnidx,:),1));
    excludedata.geneid(i) = join(genelist);
end

save('excludedata.mat','excludedata');