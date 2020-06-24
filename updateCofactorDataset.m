%% updatetmpDataset

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
for i = 1:length(CofactorUniProt.cofactor)
    tmpcofactor = CofactorUniProt.cofactor(i);
    if ismember(tmpcofactor,uniprotold)
        tmpcofactor = uniprotnew(ismember(uniprotold,tmpcofactor));
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

save('CofactorDataset.mat','CofactorDataset');




