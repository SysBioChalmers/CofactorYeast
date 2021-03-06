%% updatetmpDataset
function updateCofactorDataset
% update the collected pdbe cofactor dataset by adding more information
% from other databases, e.g., UniProt and BRENDA.

load('pdbeCofactor.mat');
load('CofactorUniProt.mat');

[~,txt,~] = xlsread('cofactorID.xlsx','PDBe');
pdbeold = txt(2:end,1);
pdbenew = txt(2:end,2);
clear txt;
[~,txt,~] = xlsread('cofactorID.xlsx','UniProt');
uniprotold = txt(2:end,1);
uniprotnew = txt(2:end,2);
clear txt;

% Convert pdbeCofactor to tmpDataset.
tmpDataset = struct;
tmpDataset.cofactor = cell(0,1);
tmpDataset.protein = cell(0,1);
tmpDataset.copy = zeros(0,1);
tmpDataset.source = cell(0,1);

for i = 1:length(pdbeCofactor.type)
    tmp_cofactor = split(pdbeCofactor.type(i),'; ');
    tmp_copy = split(pdbeCofactor.copy(i),'; ');
    tmp_copy = cellfun(@(x) str2double(x),tmp_copy,'UniformOutput',false);
    tmp_copy = cell2mat(tmp_copy);
    tmp_copy = round(tmp_copy,4);
    tmp_source = pdbeCofactor.note{i};
    if strcmp(tmp_source(1:6),'PDBid:')
        tmp_source = tmp_source(7:end);
        tmp_source = strsplit(tmp_source,'; ')';
        tmp_source = cellfun(@(x) strcat('PDBid:',x),tmp_source,'UniformOutput',false);
    else
        tmp_source = repelem({tmp_source},length(tmp_cofactor))';
    end
    tmpDataset.protein = [tmpDataset.protein;repelem(pdbeCofactor.protein(i),length(tmp_cofactor))'];
    tmpDataset.cofactor = [tmpDataset.cofactor;tmp_cofactor];
    tmpDataset.copy = [tmpDataset.copy;tmp_copy];
    tmpDataset.source = [tmpDataset.source;tmp_source];
end

% Change cofactor id and remain the cofactors of interest
CofactorDataset = struct;
CofactorDataset.cofactor = cell(0,1);
CofactorDataset.protein = cell(0,1);
CofactorDataset.copy = zeros(0,1);
CofactorDataset.source = cell(0,1);
for i = 1:length(tmpDataset.cofactor)
    tmpcofactor = tmpDataset.cofactor(i);
    if ismember(tmpcofactor,pdbeold)
        CofactorDataset.cofactor = [CofactorDataset.cofactor;pdbenew(ismember(pdbeold,tmpcofactor))];
        CofactorDataset.protein = [CofactorDataset.protein;tmpDataset.protein(i)];
        CofactorDataset.copy = [CofactorDataset.copy;tmpDataset.copy(i)];
        CofactorDataset.source = [CofactorDataset.source;tmpDataset.source(i)];
    end
end

% add new rows into CofactorDataset from UniProt database
% but iron and copper will be added in another manual file
cufe_list = {'ISC_2FE2S' 'ISC_3FE4S' 'ISC_4FE4S' 'FE_III' 'FE_II' 'HEME_A' 'HEME_B' 'HEME_C' 'CU_I' 'CU_II'};

for i = 1:length(CofactorUniProt.cofactor)
    tmpcofactor = CofactorUniProt.cofactor(i);
    if ismember(tmpcofactor,uniprotold)
        tmpcofactor = uniprotnew(ismember(uniprotold,tmpcofactor));
        if ~ismember(tmpcofactor,cufe_list)
            tmpprotein = CofactorUniProt.protein(i);
            if ~any(ismember(CofactorDataset.cofactor,tmpcofactor)&ismember(CofactorDataset.protein,tmpprotein))
            % if the cofactor of the protein has been reported in pdbe then it
            % would not be added.
                CofactorDataset.cofactor = [CofactorDataset.cofactor;tmpcofactor];
                CofactorDataset.protein = [CofactorDataset.protein;tmpprotein];
                CofactorDataset.copy = [CofactorDataset.copy;CofactorUniProt.copy(i)];
                CofactorDataset.source = [CofactorDataset.source;CofactorUniProt.source(i)];
            end
        end
    end
end


%% add manual_cofactor
% BRENDA
[num,txt,~] = xlsread('manual_cofactor.xlsx','BRENDA');

manualProtein = txt(2:end,1);
manualCofactor = txt(2:end,2);
manualSource = txt(2:end,4);
manualCopy = num;
clear num txt;
for i = 1:length(manualProtein)
    tmpprotein = manualProtein(i);
    tmpcofactor = manualCofactor(i);
    idx = ismember(CofactorDataset.cofactor,tmpcofactor)&ismember(CofactorDataset.protein,tmpprotein);
    if ~any(idx)
        CofactorDataset.cofactor = [CofactorDataset.cofactor;tmpcofactor];
        CofactorDataset.protein = [CofactorDataset.protein;tmpprotein];
        CofactorDataset.copy = [CofactorDataset.copy;manualCopy(i)];
        CofactorDataset.source = [CofactorDataset.source;manualSource(i)];
    else% if the cofactor of the protein has been collected 
        % then the copy number will be updated if the new one is larger.
        if manualCopy(i) > CofactorDataset.copy(idx)
            CofactorDataset.copy(idx) = manualCopy(i);
            CofactorDataset.source(idx) = manualSource(i);
        elseif manualCopy(i) == CofactorDataset.copy(idx)
            CofactorDataset.copy(idx) = manualCopy(i);
            CofactorDataset.source(idx) = {strcat(CofactorDataset.source{idx},';',manualSource{i})};
        end
    end
end

% Zn
[num,txt,~] = xlsread('manual_cofactor.xlsx','ZN');

manualProtein = txt(2:end,1);
manualCofactor = txt(2:end,2);
manualSource = txt(2:end,4);
manualCopy = num;
clear num txt;
for i = 1:length(manualProtein)
    tmpprotein = manualProtein(i);
    tmpcofactor = manualCofactor(i);
    idx = ismember(CofactorDataset.cofactor,tmpcofactor)&ismember(CofactorDataset.protein,tmpprotein);
    if ~any(idx)
        CofactorDataset.cofactor = [CofactorDataset.cofactor;tmpcofactor];
        CofactorDataset.protein = [CofactorDataset.protein;tmpprotein];
        CofactorDataset.copy = [CofactorDataset.copy;manualCopy(i)];
        CofactorDataset.source = [CofactorDataset.source;manualSource(i)];
    else% if the cofactor of the protein has been collected 
        % then the copy number will be updated if the new one is larger.
        if manualCopy(i) > CofactorDataset.copy(idx)
            CofactorDataset.copy(idx) = manualCopy(i);
            CofactorDataset.source(idx) = manualSource(i);
        elseif manualCopy(i) == CofactorDataset.copy(idx)
            CofactorDataset.copy(idx) = manualCopy(i);
            CofactorDataset.source(idx) = {strcat(CofactorDataset.source{idx},';',manualSource{i})};
        end
    end
end

% CU
[num,txt,~] = xlsread('manual_cofactor.xlsx','CU');

manualProtein = txt(2:end,1);
manualCofactor = txt(2:end,2);
manualSource = txt(2:end,4);
manualCopy = num;
clear num txt;
for i = 1:length(manualProtein)
    tmpprotein = manualProtein(i);
    tmpcofactor = manualCofactor(i);
    idx = ismember(CofactorDataset.cofactor,tmpcofactor)&ismember(CofactorDataset.protein,tmpprotein);
    if ~any(idx)
        CofactorDataset.cofactor = [CofactorDataset.cofactor;tmpcofactor];
        CofactorDataset.protein = [CofactorDataset.protein;tmpprotein];
        CofactorDataset.copy = [CofactorDataset.copy;manualCopy(i)];
        CofactorDataset.source = [CofactorDataset.source;manualSource(i)];
    else% if the cofactor of the protein has been collected 
        % then the copy number will be updated if the new one is larger.
        if manualCopy(i) > CofactorDataset.copy(idx)
            CofactorDataset.copy(idx) = manualCopy(i);
            CofactorDataset.source(idx) = manualSource(i);
        elseif manualCopy(i) == CofactorDataset.copy(idx)
            CofactorDataset.copy(idx) = manualCopy(i);
            CofactorDataset.source(idx) = {strcat(CofactorDataset.source{idx},';',manualSource{i})};
        end
    end
end

% FE
[num,txt,~] = xlsread('manual_cofactor.xlsx','FE');

manualProtein = txt(2:end,1);
manualCofactor = txt(2:end,2);
manualSource = txt(2:end,4);
manualCopy = num;
clear num txt;
for i = 1:length(manualProtein)
    tmpprotein = manualProtein(i);
    tmpcofactor = manualCofactor(i);
    idx = ismember(CofactorDataset.cofactor,tmpcofactor)&ismember(CofactorDataset.protein,tmpprotein);
    if ~any(idx)
        CofactorDataset.cofactor = [CofactorDataset.cofactor;tmpcofactor];
        CofactorDataset.protein = [CofactorDataset.protein;tmpprotein];
        CofactorDataset.copy = [CofactorDataset.copy;manualCopy(i)];
        CofactorDataset.source = [CofactorDataset.source;manualSource(i)];
    else% if the cofactor of the protein has been collected 
        % then the copy number will be updated if the new one is larger.
        if manualCopy(i) > CofactorDataset.copy(idx)
            CofactorDataset.copy(idx) = manualCopy(i);
            CofactorDataset.source(idx) = manualSource(i);
        elseif manualCopy(i) == CofactorDataset.copy(idx)
            CofactorDataset.copy(idx) = manualCopy(i);
            CofactorDataset.source(idx) = {strcat(CofactorDataset.source{idx},';',manualSource{i})};
        end
    end
end

% Exclude cofactors that involved in the synthesise of itself
% 1. based on literatures
[~,txt,~] = xlsread('manual_cofactor.xlsx','Exclude');

manualProtein = txt(2:end,1);
manualCofactor = txt(2:end,2);
manualNote = txt(2:end,3);
clear txt;
for i = 1:length(manualProtein)
    tmpprotein = manualProtein(i);
    tmpcofactor = manualCofactor(i);
    
    idx = ismember(CofactorDataset.cofactor,tmpcofactor)&ismember(CofactorDataset.protein,tmpprotein);
    CofactorDataset.copy(idx) = 0; % set to zero.
    CofactorDataset.source(idx) = manualNote(i);
end
% 2. based on yeast8
load('excludedata.mat');
for i = 1:length(excludedata.cofactorid)
    tmpcofactor = excludedata.cofactorid(i);
    tmpproteinlist = split(excludedata.geneid(i));
    for j = 1:length(tmpproteinlist)
        tmpprotein = tmpproteinlist(j);
        idx = ismember(CofactorDataset.cofactor,tmpcofactor)&ismember(CofactorDataset.protein,tmpprotein);
        CofactorDataset.copy(idx) = 0; % set to zero.
    end
end

idx = CofactorDataset.copy ~= 0;
CofactorDataset.cofactor = CofactorDataset.cofactor(idx);
CofactorDataset.copy = CofactorDataset.copy(idx);
CofactorDataset.protein = CofactorDataset.protein(idx);
CofactorDataset.source = CofactorDataset.source(idx);


%% Manual adjust

% delete FE_II of Ole1 collected from uniprot
idx = ismember(CofactorDataset.cofactor,'FE_II') & ismember(CofactorDataset.protein,'YGL055W') & contains(CofactorDataset.source,'UniProt');
CofactorDataset.cofactor = CofactorDataset.cofactor(~idx);
CofactorDataset.copy = CofactorDataset.copy(~idx);
CofactorDataset.protein = CofactorDataset.protein(~idx);
CofactorDataset.source = CofactorDataset.source(~idx);

% delete ISC_4FE4S of Ilv3 as it was more likely to be ISC_2FE2S (PMID: 21987576)
idx = ismember(CofactorDataset.cofactor,'ISC_4FE4S') & ismember(CofactorDataset.protein,'YJR016C');
CofactorDataset.cofactor = CofactorDataset.cofactor(~idx);
CofactorDataset.copy = CofactorDataset.copy(~idx);
CofactorDataset.protein = CofactorDataset.protein(~idx);
CofactorDataset.source = CofactorDataset.source(~idx);
if ~any(ismember(CofactorDataset.cofactor,'ISC_2FE2S') & ismember(CofactorDataset.protein,'YJR016C'))
    CofactorDataset.cofactor = [CofactorDataset.cofactor;{'ISC_2FE2S'}];
    CofactorDataset.protein = [CofactorDataset.protein;{'YJR016C'}];
    CofactorDataset.copy = [CofactorDataset.copy;1];
    CofactorDataset.source = [CofactorDataset.source;{'PMID: 21987576'}];
end

%% Add information for homologue proteins
% assume that both homologue proteins have the same cofactors and copy
% numbers.
% load yeast homologue genes
[~,txt,~] = xlsread('Yeast_homologue_genes.xlsx');%obtained 2020.01.02
list_gene = txt(2:end,2);
list_homologuegene = txt(2:end,5);
clear txt;
for i = 1:length(list_gene)
    gene_1 = list_gene(i);
    gene_2 = list_homologuegene(i);
    if ismember(gene_1,CofactorDataset.protein) || ismember(gene_2,CofactorDataset.protein)
        idx = [find(ismember(CofactorDataset.protein,gene_1));find(ismember(CofactorDataset.protein,gene_2))];
        cf = CofactorDataset.cofactor(idx);
        cp = CofactorDataset.copy(idx);
        unq_cf = unique(cf);
        for j = 1:length(unq_cf)
            max_cp = max(cp(ismember(cf,unq_cf(j))));
            idx_1 = ismember(CofactorDataset.protein,gene_1) & ismember(CofactorDataset.cofactor,unq_cf(j));
            if ~any(idx_1)
                CofactorDataset.protein = [CofactorDataset.protein;gene_1];
                CofactorDataset.cofactor = [CofactorDataset.cofactor;unq_cf(j)];
                CofactorDataset.copy = [CofactorDataset.copy;max_cp];
                CofactorDataset.source = [CofactorDataset.source;strcat('Paralog:',gene_2)];
            else
                cp_tmp = CofactorDataset.copy(idx_1);
                if cp_tmp < max_cp
                    CofactorDataset.copy(idx_1) = max_cp;
                    CofactorDataset.source(idx_1) = strcat(CofactorDataset.source(idx_1),';Copy adopted from paralog:',gene_2);
                end
            end
            
            idx_2 = ismember(CofactorDataset.protein,gene_2) & ismember(CofactorDataset.cofactor,unq_cf(j));
            if ~any(idx_2)
                CofactorDataset.protein = [CofactorDataset.protein;gene_2];
                CofactorDataset.cofactor = [CofactorDataset.cofactor;unq_cf(j)];
                CofactorDataset.copy = [CofactorDataset.copy;max_cp];
                CofactorDataset.source = [CofactorDataset.source;strcat('Paralog:',gene_1)];
            else
                cp_tmp = CofactorDataset.copy(idx_2);
                if cp_tmp < max_cp
                    CofactorDataset.copy(idx_2) = max_cp;
                    CofactorDataset.source(idx_2) = strcat(CofactorDataset.source(idx_2),';Copy adopted from paralog:',gene_1);
                end
            end
        end
    end
end


save('CofactorDataset.mat','CofactorDataset');




