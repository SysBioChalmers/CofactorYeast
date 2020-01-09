%% updateCofactorDataset

% update the collected pdbe cofactor dataset by adding more information
% from other databases, e.g., UniProt and BRENDA.

load('pdbeCofactor.mat');

Cofactor = pdbeCofactor;

% add new rows into Cofactor.










% Convert Cofactor to CofactorDataset.
CofactorDataset = struct;
CofactorDataset.cofactor = cell(0,1);
CofactorDataset.protein = cell(0,1);
CofactorDataset.copy = zeros(0,1);
CofactorDataset.source = cell(0,1);

for i = 1:length(Cofactor.type)
    tmp_cofactor = split(Cofactor.type(i),'; ');
    tmp_copy = split(Cofactor.copy(i),'; ');
    tmp_copy = cellfun(@(x) str2double(x),tmp_copy,'UniformOutput',false);
    tmp_copy = cell2mat(tmp_copy);
    tmp_copy = round(tmp_copy,4);
    tmp_source = Cofactor.note{i};
    if strcmp(tmp_source(1:6),'PDBid:')
        tmp_source = tmp_source(7:end);
        tmp_source = strsplit(tmp_source,'; ')';
        tmp_source = cellfun(@(x) strcat('PDBid:',x),tmp_source,'UniformOutput',false);
    else
        tmp_source = repelem({tmp_source},length(tmp_cofactor))';
    end
    CofactorDataset.protein = [CofactorDataset.protein;repelem(Cofactor.protein(i),length(tmp_cofactor))'];
    CofactorDataset.cofactor = [CofactorDataset.cofactor;tmp_cofactor];
    CofactorDataset.copy = [CofactorDataset.copy;tmp_copy];
    CofactorDataset.source = [CofactorDataset.source;tmp_source];
end

save('CofactorDataset.mat','CofactorDataset');