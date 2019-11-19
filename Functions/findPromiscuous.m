%% Search for promiscuous enzymes in a given GEM
%  If a protein can catalyze 2 or more different reactions, then it is
%  defined as the promiscuous protein. Besides, only the protein that
%  participates as a subunit in enzyme complex will be collected. That is, 
%  if a protein itself can catalyze multiple reactions then it would not be
%  collected here.
function promiscuous = findPromiscuous(model)

promiscuous = struct();
promiscuous.protein = cell(0,1);
promiscuous.reaction = cell(0,1);

for i = 1:length(model.genes)
    genename = model.genes(i);
    geneid = strcat('x(',num2str(i),')');
    rxn_idx = find(contains(model.rules,geneid));
    
    if length(rxn_idx) > 1
        rules_list = model.rules(rxn_idx);
        rules = strjoin(rules_list');
        if contains(rules,'&')
            rxn_list = model.rxns(rxn_idx);
            if length(rxn_idx) > 2
                rxns = strjoin(rxn_list');
                promiscuous.protein = [promiscuous.protein; genename];
                promiscuous.reaction = [promiscuous.reaction; {rxns}];
            else
                rxn_tmp = cellfun(@(x) strrep(x,'_fwd',''),rxn_list,'UniformOutput',false);
                rxn_tmp = cellfun(@(x) strrep(x,'_rvs',''),rxn_tmp,'UniformOutput',false);
                if ~strcmp(rxn_tmp{1},rxn_tmp{2})
                    rxns = strjoin(rxn_list');
                    promiscuous.protein = [promiscuous.protein; genename];
                    promiscuous.reaction = [promiscuous.reaction; {rxns}];
                end
            end
        end
    end
end


