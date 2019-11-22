%% addTranslationRxns 
function model = addComplexDilutionRxns(model)

% find all complexes
idx = contains(model.mets,'_enzyme');
cmplx_list = model.mets(idx);

for i = 1:length(cmplx_list)
    disp(['Adding complex dilution:' num2str(i) '/' num2str(length(cmplx_list))]);
    
    cmplxid = cmplx_list(i);
    rxnid = strcat(cell2mat(cmplxid),'_dilution');
    
    metlist = cmplxid;
    coeflist = -1;
    model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);
end